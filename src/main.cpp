#include <iostream>
#include <fstream>
#include <nanogui/nanogui.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#include "misc/getopt.h" // getopt for windows
#else
#include <getopt.h>
#include <unistd.h>
#endif
#include <unordered_set>
#include <stdlib.h> // atoi for getopt inputs

#include "CGL/CGL.h"
#include "collision/plane.h"
#include "collision/sphere.h"
#include "cloth.h"
#include "clothSimulator.h"
#include "json.hpp"
#include "misc/file_utils.h"

typedef uint32_t gid_t;

using namespace std;
using namespace nanogui;

using json = nlohmann::json;

#define msg(s) cerr << "[ClothSim] " << s << endl;

const string SPHERE = "sphere";
const string PLANE = "plane";
const string CLOTH = "cloth";
const string PARTICLES = "particles";
const string BOUNDARY = "bounds";

const unordered_set<string> VALID_KEYS = {SPHERE, PLANE, CLOTH, PARTICLES, BOUNDARY};

ClothSimulator *app = nullptr;
GLFWwindow *window = nullptr;
Screen *screen = nullptr;

unsigned int fbo;
unsigned int rectVAO;
unsigned int textureColorbuffer;
GLShader pixelate_shader;

Vector3D bgColor;

// from https://www.youtube.com/watch?v=QQ3jr-9Rc1o
float fullRect[] = {
    // Coords    // texCoords
     1.0f, -1.0f,  1.0f, 0.0f,
    -1.0f, -1.0f,  0.0f, 0.0f,
    -1.0f,  1.0f,  0.0f, 1.0f,

     1.0f,  1.0f,  1.0f, 1.0f,
     1.0f, -1.0f,  1.0f, 0.0f,
    -1.0f,  1.0f,  0.0f, 1.0f
};

void error_callback(int error, const char* description) {
  puts(description);
}

void onResize(int width, int height) {
    // generate texture
//    width = width * 3 / 4;
    glGenTextures(1, &textureColorbuffer);
    glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, height, width, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // Prevents edge bleeding
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Prevents edge bleeding
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);

    // attach it to currently bound framebuffer object
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);
    
    unsigned int rbo;
    glGenRenderbuffers(1, &rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, height, width);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
}

void createGLContexts(string project_root) {
  if (!glfwInit()) {
    return;
  }

  glfwSetTime(0);

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  glfwWindowHint(GLFW_SAMPLES, 0);
  glfwWindowHint(GLFW_RED_BITS, 8);
  glfwWindowHint(GLFW_GREEN_BITS, 8);
  glfwWindowHint(GLFW_BLUE_BITS, 8);
  glfwWindowHint(GLFW_ALPHA_BITS, 8);
  glfwWindowHint(GLFW_STENCIL_BITS, 8);
  glfwWindowHint(GLFW_DEPTH_BITS, 24);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create a GLFWwindow object
  window = glfwCreateWindow(800, 800, "Cloth Simulator", nullptr, nullptr);
  if (window == nullptr) {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return;
  }
  glfwMakeContextCurrent(window);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    throw std::runtime_error("Could not initialize GLAD!");
  }
  glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM

  glClearColor(0.2f, 0.25f, 0.3f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // Create a nanogui screen and pass the glfw pointer to initialize
  screen = new Screen();
  screen->initialize(window, true);

  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  glViewport(0, 0, width, height);
  glfwSwapInterval(1);
  glfwSwapBuffers(window);
    
    // post processing initialization
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    onResize(width, height);
    
    unsigned int rectVBO;
    glGenVertexArrays(1, &rectVAO);
    glGenBuffers(1, &rectVBO);
    glBindVertexArray(rectVAO);
    glBindBuffer(GL_ARRAY_BUFFER, rectVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(fullRect), &fullRect, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    
    pixelate_shader.initFromFiles("pixelate", project_root + "/shaders/pixelate.vert",
        project_root + "/shaders/pixelate.frag");
    
    pixelate_shader.setUniform("screenTexture", 0);
}

void setGLFWCallbacks() {
  glfwSetCursorPosCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->cursorPosCallbackEvent(x, y)) {
      app->cursorPosCallbackEvent(x / screen->pixelRatio(),
                                  y / screen->pixelRatio());
    }
  });

  glfwSetMouseButtonCallback(
      window, [](GLFWwindow *, int button, int action, int modifiers) {
        if (!screen->mouseButtonCallbackEvent(button, action, modifiers) ||
            action == GLFW_RELEASE) {
          app->mouseButtonCallbackEvent(button, action, modifiers);
        }
      });

  glfwSetKeyCallback(
      window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
        if (!screen->keyCallbackEvent(key, scancode, action, mods)) {
          app->keyCallbackEvent(key, scancode, action, mods);
        }
      });

  glfwSetCharCallback(window, [](GLFWwindow *, unsigned int codepoint) {
    screen->charCallbackEvent(codepoint);
  });

  glfwSetDropCallback(window,
                      [](GLFWwindow *, int count, const char **filenames) {
                        screen->dropCallbackEvent(count, filenames);
                        app->dropCallbackEvent(count, filenames);
                      });

  glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->scrollCallbackEvent(x, y)) {
      app->scrollCallbackEvent(x, y);
    }
  });

  glfwSetFramebufferSizeCallback(window,
                                 [](GLFWwindow *, int width, int height) {
                                   screen->resizeCallbackEvent(width, height);
                                   app->resizeCallbackEvent(width, height);
                                  onResize(width, height);
                                 });
}

void usageError(const char *binaryName) {
  printf("Usage: %s [options]\n", binaryName);
  printf("Required program options:\n");
  printf("  -f     <STRING>    Filename of scene\n");
  printf("  -r     <STRING>    Project root.\n");
  printf("                     Should contain \"shaders/Default.vert\".\n");
  printf("                     Automatically searched for by default.\n");
  printf("  -a     <INT>       Sphere vertices latitude direction.\n");
  printf("  -o     <INT>       Sphere vertices longitude direction.\n");
  printf("\n");
  exit(-1);
}

void incompleteObjectError(const char *object, const char *attribute) {
  cout << "Incomplete " << object << " definition, missing " << attribute << endl;
  exit(-1);
}

bool loadObjectsFromFile(string filename, Cloth *cloth, ClothParameters *cp, vector<CollisionObject *>* objects, int sphere_num_lat, int sphere_num_lon) {
  // Read JSON from file
  ifstream i(filename);
  if (!i.good()) {
    return false;
  }
  json j;
  i >> j;

  std::cout << "loaded json" << "\n";
    
    auto boundary = j.find(BOUNDARY);
    if (boundary == j.end()) {
        incompleteObjectError("JSON", "boundary data");
    }

  // Loop over objects in scene
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    string key = it.key();

    // Check that object is valid
    unordered_set<string>::const_iterator query = VALID_KEYS.find(key);
    if (query == VALID_KEYS.end()) {
      cout << "Invalid scene object found: " << key << endl;
      exit(-1);
    }

    // Retrieve object
    json object = it.value();

    // Parse object depending on type (cloth, sphere, or plane)
    if (key == PARTICLES) {
        //TODO: how to do this better?
        cloth->frames_per_sec = 90;
        cloth->simulation_steps = 30;
        
      for (auto& particleData : object) {
          int particleCount;
          Vector3D spawnPos;
          Vector3D spawnExtents;
          ParticleProperties properties;
          
          auto temp = particleData.find("count");
          if (temp != particleData.end()) {
              particleCount = *temp;
          } else {
              incompleteObjectError("particles", "count");
          }
          
          temp = particleData.find("spawnPos");
          if (temp != particleData.end()) {
              vector<double> vec_pos = *temp;
              spawnPos = Vector3D(vec_pos[0], vec_pos[1], vec_pos[2]);
          } else {
              incompleteObjectError("particles", "spawnPos");
          }
          
          temp = particleData.find("spawnExtents");
          if (temp != particleData.end()) {
              vector<double> vec_extents = *temp;
              spawnExtents = Vector3D(vec_extents[0], vec_extents[1], vec_extents[2]);
          } else {
              incompleteObjectError("particles", "spawnExtents");
          }
          
          temp = particleData.find("particleMass");
          if (temp != particleData.end()) {
              properties.mass = *temp;
          } else {
              incompleteObjectError("particles", "particleMass");
          }
          
          temp = particleData.find("particleRadius");
          if (temp != particleData.end()) {
              properties.radius = *temp;
          } else {
              incompleteObjectError("particles", "particleRadius");
          }
          temp = particleData.find("collRadius");
          if (temp != particleData.end()) {
              properties.collRadius = *temp;
          }
          else {
              incompleteObjectError("particles", "particleRadius");
          }

          temp = particleData.find("externalForces");
          if (temp != particleData.end()) {
              vector<double> vec_forces = *temp;
              properties.external_forces = Vector3D(vec_forces[0], vec_forces[1], vec_forces[2]);
          } else {
              incompleteObjectError("particles", "externalForces");
          }
          
          temp = particleData.find("velocity");
          properties.velocity = Vector3D(0);
          if (temp != particleData.end()) {
              vector<double> vec_vel = *temp;
              properties.velocity = Vector3D(vec_vel[0], vec_vel[1], vec_vel[2]);
          }
          
          temp = particleData.find("velocityColor");
          properties.velocityColor = Vector3D(0);
          if (temp != particleData.end()) {
              vector<double> vec_vel = *temp;
              properties.velocityColor = Vector3D(vec_vel[0], vec_vel[1], vec_vel[2]);
          }
          
          temp = particleData.find("pinned");
          properties.pinned = false;
          if (temp != particleData.end() && *(temp) == 1) {
              properties.pinned = true;
          }
          
          temp = particleData.find("partColl");
          if (temp != particleData.end()) {
              if (*(temp) == 1) {
                  properties.particle_collisions = true;
              }
              else {
                  properties.particle_collisions = false;
              }
          }
          else {
              incompleteObjectError("particles", "partColl");
          }
          
          temp = particleData.find("shaderType");
          if (temp != particleData.end()) {
              properties.shaderType = *temp;
          }else {
              incompleteObjectError("particles", "shaderType");
          }
          
          temp = particleData.find("particleAveragingFactor");
          properties.particleAveragingFactor = 0;
          if (temp != particleData.end()) {
              properties.particleAveragingFactor = *temp;
          }
          
          temp = particleData.find("particleAveragingDist");
          properties.particleAveragingDist = 0;
          if (temp != particleData.end()) {
              properties.particleAveragingDist = *temp;
          }
          
          temp = particleData.find("particleAveragingBrightness");
          properties.particleAveragingBrightness = 1;
          if (temp != particleData.end()) {
              properties.particleAveragingBrightness = *temp;
          }
          
          temp = particleData.find("particleColor");
          if (temp != particleData.end()) {
              vector<double> vec_color = *temp;
              properties.color = Vector3D(vec_color[0], vec_color[1], vec_color[2]);
          } else {
              incompleteObjectError("particles", "particleColor");
          }

          temp = particleData.find("transformations");
          if (temp != particleData.end()) {
              vector<double> vec_transform = *temp;
              for (int val : vec_transform) {
                  properties.collision_transformations.push_back(val);
              }
          }
          else {
              incompleteObjectError("particles", "transformations");
          }
          
          forceLaw curLaw;
          temp = particleData.find("receivedForces");
          if (temp != particleData.end()) {
              for (auto& forceData : particleData["receivedForces"]) {
                  std::cout << "start forces" << "\n";
                  
                  temp = forceData.find("law");
                  
                  if (temp != forceData.end()) {
                      auto lawName = *temp;
                      
                      if (lawName == "r2") {
                          curLaw = &r2_law;
                      } else if (lawName == "r4") {
                          curLaw = &r4_law;
                      } else if (lawName == "cross") {
                          curLaw = &cross_law;
                      }else if (lawName == "fireForce") {
                          curLaw = &fire_force;
                      }else {
                          incompleteObjectError("force", "law1");
                      }
                                        
                      properties.force_laws.emplace_back(curLaw);
                  } else {
                      incompleteObjectError("force", "law2");
                  }
                  
                  std::cout << "end forces" << "\n";
                  
                  temp = forceData.find("strengths");
                  if (temp != forceData.end()) {
                      vector<float> vec_strengths = *temp;
                      properties.strengths.push_back(vec_strengths);
                  } else {
                      incompleteObjectError("force", "strengths");
                  }

                  temp = forceData.find("localized");
                  if (temp != forceData.end()) {
                      if (*temp == 1) {
                          properties.localized.push_back(true);
                      }
                      else {
                          properties.localized.push_back(false);
                      }
                      
                  }
                  else {
                      incompleteObjectError("force", "strengths");
                  }
              }
              
              cloth->spawnParticles(particleCount, spawnPos, spawnExtents, properties);
          } else {
              incompleteObjectError("particles", "receivedForces");
          }
      }

      std::cout << "set up particles" << "\n";
        
    } else if (key == BOUNDARY) {
        double marchingSize;
        double physicsBuffer;
        int sideCellCounts[3];
        
        auto temp = object.find("marchingSize");
        if (temp != object.end()) {
            marchingSize = *temp;
        } else {
            incompleteObjectError("bounds", "marchingSize");
        }
        
        temp = object.find("marchingCount");
        if (temp != object.end()) {
            vector<int> numCells = *temp;
            sideCellCounts[0] = numCells[0];
            sideCellCounts[1] = numCells[1];
            sideCellCounts[2] = numCells[2];
        } else {
            incompleteObjectError("bounds", "marchingCount");
        }
        
        temp = object.find("noTop");
        bool noTop = false;
        if (temp != object.end() && *temp == 1) {
            noTop = true;
        }
        
        temp = object.find("damping");
        double damping = 0.01;
        if (temp != object.end()) {
            damping = *temp;
        }
        cloth->damping = damping;
        
        temp = object.find("physicsBuffer");
        if (temp != object.end()) {
            physicsBuffer = *temp;
        } else {
            incompleteObjectError("bounds", "physicsBuffer");
        }
        
        temp = object.find("inSphere");
        cloth->isInSphere = false;
        if (temp != object.end() && *temp == 1) {
            cloth->isInSphere = true;
        }
        
        temp = object.find("useDensity");
        cloth->useDensity = false;
        if (temp != object.end() && *temp == 1) {
            cloth->useDensity = true;
        }
        
        bgColor = Vector3D(0.25);
        temp = object.find("bg");
        if (temp != object.end()) {
            vector<double> vec_bg = *temp;
            bgColor = Vector3D(vec_bg[0], vec_bg[1], vec_bg[2]);
        }
        
        cloth->initMarchingCubes(sideCellCounts[0], sideCellCounts[1], sideCellCounts[2], marchingSize, physicsBuffer, noTop);
    } else if (key == SPHERE) {
      Vector3D origin;
      double radius, friction;

      auto it_origin = object.find("origin");
      if (it_origin != object.end()) {
        vector<double> vec_origin = *it_origin;
        origin = Vector3D(vec_origin[0], vec_origin[1], vec_origin[2]);
      } else {
        incompleteObjectError("sphere", "origin");
      }

      auto it_radius = object.find("radius");
      if (it_radius != object.end()) {
        radius = *it_radius;
      } else {
        incompleteObjectError("sphere", "radius");
      }

      auto it_friction = object.find("friction");
      if (it_friction != object.end()) {
        friction = *it_friction;
      } else {
        incompleteObjectError("sphere", "friction");
      }

      for (int k = 0; k < 3; k++) {
          float a = (rand() % 100) / 100.0 - 0.5;
          float b = 2 * (rand() % 100) / 100.0 - 1.4;
          float c = (rand() % 100) / 100.0 - 0.5;
          Sphere* s = new Sphere(origin + Vector3D(2 * a, 4 * b, 2 * c), radius, friction, sphere_num_lat, sphere_num_lon);
          objects->push_back(s);
      }
        
      
    } else if (key == PLANE) { // PLANE
      Vector3D point, normal;
      double friction;

      auto it_point = object.find("point");
      if (it_point != object.end()) {
        vector<double> vec_point = *it_point;
        point = Vector3D(vec_point[0], vec_point[1], vec_point[2]);
      } else {
        incompleteObjectError("plane", "point");
      }

      auto it_normal = object.find("normal");
      if (it_normal != object.end()) {
        vector<double> vec_normal = *it_normal;
        normal = Vector3D(vec_normal[0], vec_normal[1], vec_normal[2]);
      } else {
        incompleteObjectError("plane", "normal");
      }

      auto it_friction = object.find("friction");
      if (it_friction != object.end()) {
        friction = *it_friction;
      } else {
        incompleteObjectError("plane", "friction");
      }

      Plane *p = new Plane(point, normal, friction);
      objects->push_back(p);
    } else {
        cout << "Unexpected type: " << key << endl;
        exit(-1);
        
    }
  }

  i.close();
  
  return true;
}

bool is_valid_project_root(const std::string& search_path) {
    std::stringstream ss;
    ss << search_path;
    ss << "/";
    ss << "shaders/Default.vert";
    
    return FileUtils::file_exists(ss.str());
}

// Attempt to locate the project root automatically
bool find_project_root(const std::vector<std::string>& search_paths, std::string& retval) {
  
  for (std::string search_path : search_paths) {
    if (is_valid_project_root(search_path)) {
      retval = search_path;
      return true;
    }
  }
  return false;
}

int main(int argc, char **argv) {
  // Attempt to find project root
  std::vector<std::string> search_paths = {
    ".",
    "..",
    "../..",
    "../../.."
  };
  std::string project_root;
  bool found_project_root = find_project_root(search_paths, project_root);
  
  Cloth cloth;
  ClothParameters cp;
  vector<CollisionObject *> objects;
  
  int c;
  
  int sphere_num_lat = 40;
  int sphere_num_lon = 40;
  
  std::string file_to_load_from;
  bool file_specified = false;
  
  while ((c = getopt (argc, argv, "f:r:a:o:p")) != -1) {
    switch (c) {
      case 'f': {
        file_to_load_from = optarg;
        file_specified = true;
        break;
      }
      case 'r': {
        project_root = optarg;
        if (!is_valid_project_root(project_root)) {
          std::cout << "Warn: Could not find required file \"shaders/Default.vert\" in specified project root: " << project_root << std::endl;
        }
        found_project_root = true;
        break;
      }
      case 'a': {
        int arg_int = atoi(optarg);
        if (arg_int < 1) {
          arg_int = 1;
        }
        sphere_num_lat = arg_int;
        break;
      }
      case 'o': {
        int arg_int = atoi(optarg);
        if (arg_int < 1) {
          arg_int = 1;
        }
        sphere_num_lon = arg_int;
        break;
      }
      default: {
        usageError(argv[0]);
        break;
      }
    }
  }
  
  if (!found_project_root) {
    std::cout << "Error: Could not find required file \"shaders/Default.vert\" anywhere!" << std::endl;
    return -1;
  } else {
    std::cout << "Loading files starting from: " << project_root << std::endl;
  }

  if (!file_specified) { // No arguments, default initialization
    std::stringstream def_fname;
    def_fname << project_root;
    def_fname << "/scene/pinned2.json";
    file_to_load_from = def_fname.str();
  }
  bool success = loadObjectsFromFile(file_to_load_from, &cloth, &cp, &objects, sphere_num_lat, sphere_num_lon);
  if (!success) {
    std::cout << "Warn: Unable to load from file: " << file_to_load_from << std::endl;
  }

  glfwSetErrorCallback(error_callback);

  createGLContexts(project_root);

  // Initialize the ClothSimulator object
  app = new ClothSimulator(project_root, screen);
  app->loadCloth(&cloth);
  app->loadClothParameters(&cp);
  app->loadCollisionObjects(&objects);
  app->init();

  // Call this after all the widgets have been defined

  screen->setVisible(true);
  screen->performLayout();

  // Attach callbacks to the GLFW window

  setGLFWCallbacks();
    
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    onResize(width, height);
    
    app->should_pixelate = false;

  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
            
      if (app->should_pixelate) {
          // Bind the custom framebuffer
          glBindFramebuffer(GL_FRAMEBUFFER, fbo);
      } else {
          glBindFramebuffer(GL_FRAMEBUFFER, 0);
      }

    glClearColor(bgColor.x, bgColor.y, bgColor.z, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    app->drawContents();

    // Draw nanogui
    screen->drawContents();
//    screen->drawWidgets();
      
      if (app->should_pixelate) {
          // Bind the default framebuffer
          glBindFramebuffer(GL_FRAMEBUFFER, 0);
          // Draw the framebuffer rectangle
          pixelate_shader.bind();
          glBindVertexArray(rectVAO);
          glDisable(GL_DEPTH_TEST); // prevents framebuffer rectangle from being discarded
          glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
          glDrawArrays(GL_TRIANGLES, 0, 6);
      }

    glfwSwapBuffers(window);

    if (!app->isAlive()) {
      glfwSetWindowShouldClose(window, 1);
    }
  }
    
    glDeleteFramebuffers(1, &fbo);

  return 0;
}
