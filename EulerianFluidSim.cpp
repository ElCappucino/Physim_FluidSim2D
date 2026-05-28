#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL

#include <glad/glad.h>   // 1. GLAD always first
#include <GLFW/glfw3.h>  // 2. GLFW second

// 4. Standard Libraries
#include <iostream>
#include <vector>
#include <algorithm>

// 5. External Libraries (stb, glm)
#include <stb_image.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// 6. Your Custom LearnOpenGL Wrappers
#include <learnopengl/shader_m.h>

// 7. ImGui Headers
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

enum SimState
{
    TANK = 0,
    WIND_TUNNEL,
    PAINT,
    HIRES_TUNNEL
};

bool isHoldingLeftClick = false;
bool isShowingStreamline = true;

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// Parameters
float gravity = -9.81f;

int numIters = 100;
float overRelaxation = 1.9f;
float dt;
int frameNr = 0; // how many simulation step has passed
glm::vec2 ObjectaclePos;
float obstacleRadius = 0.15f;
bool pause = false;
SimState sceneNr = SimState::WIND_TUNNEL; // current state
bool showObstacle = false;
bool showStreamlines = false;
bool showPressure = false;
bool showSmoke = true;
float domainHeight, domainWidth;


class Fluid
{
public:
    float density;
    int num;
    int numCells;
    int numX, numY;
    float h;
    std::vector<float> u;
    std::vector<float> newU;
    std::vector<float> v;
    std::vector<float> newV;
    std::vector<float> p;
    std::vector<float> s;
    std::vector<float> m;
    std::vector<float> newM;



    enum FieldType
    {
        U_FIELD = 0,
        V_FIELD = 1,
        S_FIELD = 2
    };

    Fluid(float density, int numX, int numY, float h)
    {
        this->density = density;
        this->numX = numX + 2;
        this->numY = numY + 2;
        this->numCells = this->numX * this->numY;
        this->h = h;

        u.assign(numCells, 0.0f);
        v.assign(numCells, 0.0f);
        newU.assign(numCells, 0.0f);
        newV.assign(numCells, 0.0f);
        p.assign(numCells, 0.0f);
        s.assign(numCells, 0.0f);
        newM.assign(numCells, 0.0f);

        m.assign(numCells, 1.0f);
    }

    void Integrate(float dt, float gravity)
    {
        int n = this->numY;
        for (int i = 1; i < this->numY; i++)
        {
            for (int j = 1; j < this->numY - 1; j++)
            {
                if (this->s[i * n + j] != 0.0f && this->s[i * n + j - 1] != 0.0f)
                {
                    this->v[i * n + j] += gravity * dt;
                }
            }
        }
    }

    void SolveIncompressibility(int numIters, float dt)
    {
        int n = this->numY;
        float cp = this->density * this->h / dt;

        for (int iter = 0; iter < numIters; iter++)
        {
            for (int i = 1; i < this->numX - 1; i++)
            {
                for (int j = 1; j < this->numY - 1; j++)
                {
                    if (this->s[i * n + j] == 0.0f) continue;

                    float sx0 = this->s[(i - 1) * n + j];
                    float sx1 = this->s[(i + 1) * n + j];
                    float sy0 = this->s[i * n + j - 1];
                    float sy1 = this->s[i * n + j + 1];
                    float s = sx0 + sx1 + sy0 + sy1;

                    if (s == 0.0f) continue;

                    float div = this->u[(i + 1) * n + j] - this->u[i * n + j] +
                        this->v[i * n + j + 1] - this->v[i * n + j];

                    float pVar = -div / s;
                    pVar *= overRelaxation;
                    this->p[i * n + j] += cp * pVar;

                    this->u[i * n + j] -= sx0 * pVar;
                    this->u[(i + 1) * n + j] += sx1 * pVar;

                    this->v[i * n + j] -= sy0 * pVar;
                    this->v[i * n + j + 1] += sy1 * pVar;
                }
            }
        }
    }

    void Extrapolate()
    {
        int n = this->numY;
        for (int i = 0; i < this->numX; i++) {
            this->u[i * n + 0] = this->u[i * n + 1];
            this->u[i * n + n - 1] = this->u[i * n + n - 2];
        }
        for (int j = 0; j < this->numY; j++) {
            this->v[0 * n + j] = v[1 * n + j];
            this->v[(this->numX - 1) * n + j] = this->v[(this->numX - 2) * n + j];
        }
    }

    float SampleField(float x, float y, FieldType field)
    {
        int n = this->numY;
        float h1 = 1.0f / h;
        float h2 = 0.5f * h;

        x = (std::max)((std::min)(x, this->numX * this->h), this->h);
        y = (std::max)((std::min)(y, this->numY * this->h), this->h);

        float dx = 0.0f;
        float dy = 0.0f;

        std::vector<float>* f = nullptr;

        switch (field) {
        case U_FIELD: f = &this->u; dy = h2; break;
        case V_FIELD: f = &this->v; dx = h2; break;
        case S_FIELD: f = &this->m; dx = h2; dy = h2; break;
        }

        int x0 = (std::min)((int)std::floor((x - dx) * h1), this->numX - 1);
        float tx = ((x - dx) - x0 * this->h) * h1;
        int x1 = (std::min)(x0 + 1, this->numX - 1);

        int y0 = (std::min)((int)std::floor((y - dy) * h1), this->numY - 1);
        float ty = ((y - dy) - y0 * this->h) * h1;
        int y1 = (std::min)(y0 + 1, this->numY - 1);

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        float val = sx * sy * (*f)[x0 * n + y0] +
            tx * sy * (*f)[x1 * n + y0] +
            tx * ty * (*f)[x1 * n + y1] +
            sx * ty * (*f)[x0 * n + y1];

        return val;
    }

    float AvgU(int i, int j) {
        int n = this->numY;
        float u = (this->u[i * n + j - 1] + this->u[i * n + j] +
            this->u[(i + 1) * n + j - 1] + this->u[(i + 1) * n + j]) * 0.25f;
        return u;
    }

    float AvgV(int i, int j) {
        int n = this->numY;
        float u = (this->v[(i - 1) * n + j] + this->v[i * n + j] +
            this->v[(i - 1) * n + j + 1] + this->v[i * n + j + 1]) * 0.25f;
        return u;
    }

    void AdvectVel(float dt) {
        newU = this->u;
        newV = this->v;

        int n = this->numY;
        float h2 = 0.5f * this->h;

        for (int i = 1; i < this->numX; i++) {
            for (int j = 1; j < this->numY; j++) {
                // u component
                if (this->s[i * n + j] != 0.0f && this->s[(i - 1) * n + j] != 0.0f && j < this->numY - 1) {
                    float x = i * this->h;
                    float y = j * this->h + h2;
                    float currentU = this->u[i * n + j];
                    float currentV = AvgV(i, j);
                    x -= dt * currentU;
                    y -= dt * currentV;
                    this->newU[i * n + j] = SampleField(x, y, U_FIELD);
                }
                // v component
                if (this->s[i * n + j] != 0.0f && this->s[i * n + j - 1] != 0.0f && i < this->numX - 1) {
                    float x = i * this->h + h2;
                    float y = j * this->h;
                    float currentU = AvgU(i, j);
                    float currentV = this->v[i * n + j];
                    x -= dt * currentU;
                    y -= dt * currentV;
                    this->newV[i * n + j] = SampleField(x, y, V_FIELD);
                }
            }
        }
        u = this->newU;
        v = this->newV;
    }

    void AdvectSmoke(float dt) {

        this->newM = this->m;
        int n = this->numY;
        float h2 = 0.5f * h;

        for (int i = 1; i < this->numX - 1; i++) {
            for (int j = 1; j < this->numY - 1; j++) {

                if (this->s[i * n + j] != 0.0f) {
                    float currentU = (this->u[i * n + j] + this->u[(i + 1) * n + j]) * 0.5f;
                    float currentV = (this->v[i * n + j] + this->v[i * n + j + 1]) * 0.5f;
                    float x = (i * h) + h2 - (dt * currentU);
                    float y = (j * h) + h2 - (dt * currentV);

                    this->newM[i * n + j] = SampleField(x, y, S_FIELD);
                }
            }
        }
        this->m = this->newM;
    }

    void Fluid::simulate(float dt, float gravity, int numIters) {
        Integrate(dt, gravity);

        std::fill(p.begin(), p.end(), 0.0f);
        SolveIncompressibility(numIters, dt);

        Extrapolate();
        AdvectVel(dt);
        AdvectSmoke(dt);
    }
};

Fluid* fluid;

const char *vertexShaderSource = "#version 330 core\n"
    "layout (location = 0) in vec3 aPos;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
    "}\0";
const char *fragmentShaderSource = "#version 330 core\n"
    "out vec4 FragColor;\n"
    "void main()\n"
    "{\n"
    "   FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
    "}\n\0";



void setObstacle(float x, float y, bool reset) {
    float vx = 0.0f;
    float vy = 0.0f;

    // Calculate velocity of the obstacle based on mouse movement
    if (!reset) {
        vx = (x - ObjectaclePos.x) / dt;
        vy = (y - ObjectaclePos.y) / dt;
    }

    ObjectaclePos.x = x;
    ObjectaclePos.y = y;

    float r = obstacleRadius;
    int n = fluid->numY;

    // Loop through all inner cells
    for (int i = 1; i < fluid->numX - 2; i++) {
        for (int j = 1; j < fluid->numY - 2; j++) {

            // Reset cell to fluid first
            fluid->s[i * n + j] = 1.0f;

            // Calculate distance from center of cell to obstacle center
            float dx = (i + 0.5f) * fluid->h - x;
            float dy = (j + 0.5f) * fluid->h - y;

            // If the cell is inside the circle (distance squared < radius squared)
            if (dx * dx + dy * dy < r * r) {
                fluid->s[i * n + j] = 0.0f; // Mark as solid

                if (sceneNr == 2) { // Paint mode: dynamic color
                    fluid->m[i * n + j] = 0.5f + 0.5f * std::sin(0.1f * frameNr);
                }
                else {
                    fluid->m[i * n + j] = 1.0f; // Standard smoke
                }

                // Set fluid velocity to match the obstacle's velocity
                fluid->u[i * n + j] = vx;
                fluid->u[(i + 1) * n + j] = vx;
                fluid->v[i * n + j] = vy;
                fluid->v[i * n + j + 1] = vy;
            }
        }
    }
    showObstacle = true;
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // build and compile our shader program
    // ------------------------------------
    // vertex shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // link shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float vertices[] = {
        -0.5f, -0.5f, 0.0f, // left  
         0.5f, -0.5f, 0.0f, // right 
         0.0f,  0.5f, 0.0f  // top   
    }; 

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0); 

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    glBindVertexArray(0); 

    float quadVertices[] = {
        // positions        // texture coords (used for circle math)
        -1.0f,  1.0f, 0.0f,  0.0f, 1.0f,
        -1.0f, -1.0f, 0.0f,  0.0f, 0.0f,
         1.0f, -1.0f, 0.0f,  1.0f, 0.0f,

        -1.0f,  1.0f, 0.0f,  0.0f, 1.0f,
         1.0f, -1.0f, 0.0f,  1.0f, 0.0f,
         1.0f,  1.0f, 0.0f,  1.0f, 1.0f
    };

    unsigned int quadVAO, quadVBO;
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // TexCoord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Setup Parameter
    ObjectaclePos = glm::vec2(0.0f, 0.0f);
    obstacleRadius = 0.15f;
    overRelaxation = 1.9f;
    dt = 1.0f / 120.0f;
    numIters = 40;

    int res = 100;
    if (sceneNr == 0) res = 50;
    else if (sceneNr == 3) res = 200;

    domainHeight = 1.0f;
    domainWidth = (domainHeight / (float)SCR_HEIGHT) * (float)SCR_WIDTH;
    float h = domainHeight / res;

    int numX = std::floor(domainWidth / h);
    int numY = std::floor(domainHeight / h);
    float density = 1000.0f;

    if (fluid != nullptr) delete fluid;
    fluid = new Fluid(density, numX, numY, h);

    int n = fluid->numY;


    if (sceneNr == SimState::TANK) { // TANK
        for (int i = 0; i < fluid->numX; i++) {
            for (int j = 0; j < fluid->numY; j++) {
                float s = 1.0f; // fluid
                if (i == 0 || i == fluid->numX - 1 || j == 0)
                    s = 0.0f; // solid walls
                fluid->s[i * n + j] = s;
            }
        }
        gravity = -9.81f;
        showPressure = true;
        showSmoke = false;
    }
    else if (sceneNr == SimState::WIND_TUNNEL || sceneNr == SimState::HIRES_TUNNEL) { // VORTEX SHEDDING
        float inVel = 2.0f;
        for (int i = 0; i < fluid->numX; i++) {
            for (int j = 0; j < fluid->numY; j++) {
                float s = 1.0f;
                if (i == 0 || j == 0 || j == fluid->numY - 1)
                    s = 0.0f;
                fluid->s[i * n + j] = s;

                if (i == 1) fluid->u[i * n + j] = inVel;
            }
        }

        // Add a "smoke" pipe effect
        int pipeH = (int)(0.1f * fluid->numY);
        int minJ = std::floor(0.5f * fluid->numY - 0.5f * pipeH);
        int maxJ = std::floor(0.5f * fluid->numY + 0.5f * pipeH);

        for (int i = 0; i < 5; i++) { // First 5 columns
            for (int j = minJ; j < maxJ; j++) {
                fluid->m[i * n + j] = 0.0f;
            }
        }

        setObstacle(0.4f, 0.5f, true);

        gravity = 0.0f;
        showSmoke = true;

        if (sceneNr == 3) {
            dt = 1.0f / 120.0f;
            numIters = 100;
            showPressure = true;
        }
    }
    else if (sceneNr == SimState::PAINT) { // PAINT
        gravity = 0.0f;
        overRelaxation = 1.0f;
        showSmoke = true;
        obstacleRadius = 0.1f;
    }

    // Initialize streamline ui
    std::vector<float> streamlineVertices;
    unsigned int streamVAO, streamVBO;

    glGenVertexArrays(1, &streamVAO);
    glGenBuffers(1, &streamVBO);
    glBindVertexArray(streamVAO);
    glBindBuffer(GL_ARRAY_BUFFER, streamVBO);

    // Initialize with NULL; we will update this every frame
    glBufferData(GL_ARRAY_BUFFER, 100000 * sizeof(float), NULL, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    unsigned int smokeTex;
    glGenTextures(1, &smokeTex);
    glBindTexture(GL_TEXTURE_2D, smokeTex);

    // Upload the fluid->m vector to the texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, fluid->numY, fluid->numX, 0, GL_RED, GL_FLOAT, fluid->m.data());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);


    Shader circleShader("circle.vs", "circle.fs");
    Shader fluidShader("fluid.vs", "fluid.fs");

    // uncomment this call to draw in wireframe polygons.
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

        if (!pause) {
            fluid->simulate(dt, gravity, numIters);
            frameNr++;
        }

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, smokeTex);
        // Update existing texture memory with new fluid->m data
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, fluid->numY, fluid->numX, GL_RED, GL_FLOAT, fluid->m.data());

        fluidShader.use();
        fluidShader.setBool("showSmoke", showSmoke);
        fluidShader.setInt("sceneNr", sceneNr);

        
        
        // 3. Draw the full-screen quad (the background)
        glBindVertexArray(quadVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);


        // Render Streamline
        if (isShowingStreamline)
        {
            streamlineVertices.clear();
            float totalWidth = fluid->numX * fluid->h;
            float totalHeight = fluid->numY * fluid->h;

            for (int i = 1; i < fluid->numX - 1; i += 10) { // Skip every 5 cells for clarity
                for (int j = 1; j < fluid->numY - 1; j += 10) {
                    float x = (i + 0.5f) * fluid->h;
                    float y = (j + 0.5f) * fluid->h;

                    for (int n_step = 0; n_step < 15; n_step++) {
                        // Start Point
                        streamlineVertices.push_back((x / totalWidth) * 2.0f - 1.0f);
                        streamlineVertices.push_back((y / totalHeight) * 2.0f - 1.0f);
                        streamlineVertices.push_back(0.0f);

                        float u_val = fluid->SampleField(x, y, Fluid::U_FIELD);
                        float v_val = fluid->SampleField(x, y, Fluid::V_FIELD);
                        x += u_val * 0.01f;
                        y += v_val * 0.01f;

                        // End Point
                        streamlineVertices.push_back((x / totalWidth) * 2.0f - 1.0f);
                        streamlineVertices.push_back((y / totalHeight) * 2.0f - 1.0f);
                        streamlineVertices.push_back(0.0f);
                    }
                }
            }





            glUseProgram(shaderProgram); // Use a shader that outputs a solid color (e.g. Cyan)

            // Upload and Draw
            glBindBuffer(GL_ARRAY_BUFFER, streamVBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0, streamlineVertices.size() * sizeof(float), streamlineVertices.data());
            glBindVertexArray(streamVAO);
            glDrawArrays(GL_LINES, 0, streamlineVertices.size() / 3);
        }
        

        if (showObstacle)
        {
            circleShader.use();

            float scaleX = (obstacleRadius / domainWidth) * 2.0f;
            float scaleY = (obstacleRadius / domainHeight) * 2.0f;

            // Convert sim position to NDC (-1 to 1)
            float ndcX = (ObjectaclePos.x / domainWidth) * 2.0f - 1.0f;
            float ndcY = (ObjectaclePos.y / domainHeight) * 2.0f - 1.0f;

            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, glm::vec3(ndcX, ndcY, 0.0f));
            model = glm::scale(model, glm::vec3(scaleX, scaleY, 1.0f));

            // Set the uniform using the class helper
            circleShader.setMat4("model", model);
            circleShader.setVec3("obstacleColor", 0.9f, 0.9f, 0.9f);
            circleShader.setFloat("outlineThickness", 0.01f);
            circleShader.setVec3("outlineColor", 0.0f, 0.0f, 0.0f);

            glBindVertexArray(quadVAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }
        

        
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(550, 0), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(250, 300), ImGuiCond_Once);

        ImGui::Begin("Hello");
        ImGui::Text("Right click to move the obstacle");
        ImGui::Checkbox("isShowingStreamline", &isShowingStreamline);
        ImGui::Checkbox("showObstacle", &showObstacle);
        ImGui::Checkbox("Pause", &pause);
        ImGui::DragFloat("gravity", &gravity, 0.01f, -0.01f, -9.81f);
        ImGui::DragInt("Numbers of iteration", &numIters, 1, 10, 1000);
        ImGui::DragFloat("Overrelaxation", &overRelaxation, 0.01f, 0.1f, 2.0f);
        ImGui::DragFloat("Obstacle Radius", &obstacleRadius, 0.01f, 0.1f, 2.0f);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        //// draw our first triangle
        //glUseProgram(shaderProgram);
        //glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
        //glDrawArrays(GL_TRIANGLES, 0, 3);
        // glBindVertexArray(0); // no need to unbind it every time 
 
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);

    

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS)
    {
        isHoldingLeftClick = true;
    }
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_2) == GLFW_RELEASE)
    {
        isHoldingLeftClick = false;
    }

    if (isHoldingLeftClick)
    {
        // Convert screen pixels to 0.0 -> domainWidth/Height
        float simX = (mouseX / (float)SCR_WIDTH) * domainWidth;
        float simY = (1.0f - (mouseY / (float)SCR_HEIGHT)) * domainHeight; // Flip Y

        setObstacle(simX, simY, false);
    }
    
}
// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}