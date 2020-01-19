#include "../../glm/glm.hpp" // glm::vec3, glm::vec4, glm::ivec4, glm::mat4

struct triangle
{
	glm::vec3 Position[3];
	glm::vec3 Normal;
};

void computeNormal(triangle & Triangle)
{
	glm::vec3 const & a = Triangle.Position[0];
	glm::vec3 const & b = Triangle.Position[1];
	glm::vec3 const & c = Triangle.Position[2];
	Triangle.Normal = glm::normalize(glm::cross(c - a, b - a));
}

#include "../../glm/gtc/matrix_projection.hpp" // glm::perspective
#include "../../glm/gtc/matrix_transform.hpp" // glm::translate, glm::rotate, glm::scale
#include "../../glm/gtc/type_ptr.hpp" // glm::value_ptr

glm::vec4 const ClearColor = glm::vec4(glm::vec3(0.0f), 1.0f);
glm::ivec4 const Viewport = glm::ivec4(0, 0, 640, 480);

void glUniformMatrix4fv(int location, int count, bool transpose, const float * value);

int main()
{
	//triangle Triangle;
	//computeNormal(Triangle);

	float Translate = 1.0f;
	glm::vec2 Rotate;

	glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.f);
	glm::mat4 ViewTranslate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -Translate));
	glm::mat4 ViewRotateX = glm::rotate(ViewTranslate, Rotate.y, glm::vec3(-1.0f, 0.0f, 0.0f));
	glm::mat4 View = glm::rotate(ViewRotateX, Rotate.x, glm::vec3(0.0f, 1.0f, 0.0f));
	glm::mat4 Model = glm::scale(glm::mat4(1.0f), glm::vec3(0.5f));
	glm::mat4 MVP = Projection * View * Model;
	
	glUniformMatrix4fv(LocationMVP, 1, GL_FALSE, glm::value_ptr(MVP));
}
