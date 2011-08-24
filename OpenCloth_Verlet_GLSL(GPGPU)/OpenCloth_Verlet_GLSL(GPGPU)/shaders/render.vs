 
uniform int selected_index;
varying vec4 out_color;
uniform vec4 color;
void main() {
	gl_Position = ftransform(); 
	if(selected_index==gl_VertexID)
	   out_color = vec4(1)-color;
	else
	   out_color = color;
}