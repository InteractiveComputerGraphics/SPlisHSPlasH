#version 330

uniform float radius;
uniform float min_scalar;
uniform float max_scalar;
uniform mat4 projection_matrix;
uniform sampler1D samp;

in block
{
	flat vec3 mv_pos;
	flat float scalar_field;
}
In;

out vec4 out_color;
uniform vec3 color;


vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

void main(void)
{
    // calculate normal 
    vec3 n;
    n.xy = gl_PointCoord* 2.0 - vec2(1.0);    
    float mag = dot(n.xy, n.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
    n.z = sqrt(1.0-mag);

    // calculate lighting
	const vec3 light_dir = vec3(0.0, 0.0, 1.0);
    float diffuse = max(0.0, dot(light_dir, n));
 
	vec3 eye = In.mv_pos + vec3(0.0, 0.0, radius * n.z);
    vec3 halfVector = normalize( eye + light_dir);	
    float spec = pow(max(0.0, dot(n,halfVector)), 100.0);
	
	float depth = (projection_matrix[2][2] * eye.z + projection_matrix[3][2])
        / (projection_matrix[2][3] * eye.z + projection_matrix[3][3]);

    gl_FragDepth = (depth + 1.0) / 2.0;

	// modify color according to the scalar field
	vec3 hsv = rgb2hsv(color);
	float v = max(In.scalar_field-min_scalar, 0.0);
	float diff = abs(max_scalar-min_scalar);
	v = min(v/diff, 1.0);	
	vec3 fluidColor = texture(samp, v).xyz;

	// compute final color
	vec3 color_ = 0.25 * fluidColor;
    color_ += 0.7 * diffuse * fluidColor;
    color_ += 0.05 * spec * vec3(1.0);
	color_ = clamp(color_, 0.0, 1.0);
	
    out_color = vec4(color_, 1.0);
}