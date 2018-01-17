#define PARTICLESNUM 2
#define THRESHOLD 0.40
#define DT 0.1
#define STEP 0.05
#define ENTRY 0
#define EXIT 1
#define MAXSTEP 400
#define BBOX_SIZE 2.0

uniform float time;
uniform vec2 resolution;
uniform sampler2D positionSampler;
uniform sampler2D particlesSampler;
uniform sampler2D velocitySampler;
uniform vec3 eyeCoordinate;

uniform vec3 light1;
uniform vec3 light2;
uniform vec3 light3;

struct Intersection
{
    float t;  // the scalar distance
    int pIdx; // the corresponing particle's index
    int type; // entry or exit
};

struct Particle
{
    vec3 pos;
    vec3 vel;
    float radius;
    vec4 color;
    vec3 x1; // intersection 1
    vec3 x2; // intersection 2
    vec3 w;
    bool active;
};

bool sphereIntersections(vec3 c, float r, vec3 o, vec3 d, inout float t1, inout float t2)
{
    // Ray/Sphere Intersection Test
    // float A = 1;
    vec3 v = o - c;
    float B = 2.0 * dot(d, v);
    float C = dot(v, v) - r*r;

    float det = B*B - 4.0*C;

    if (det < 0.0)
        return false; // No solution

    t1 = -B - sqrt(det);
    t2 = -B + sqrt(det);

    if (t1 < 0.0 || t2 < 0.0)
        return false; // Sphere is behind the ray

    return true; 
}

bool sphereIntersections(Particle p, vec3 o, vec3 d, inout float t1, inout float t2)
{
    return sphereIntersections(p.pos, p.radius, o, d, t1, t2);
}

vec3 getCartesian(vec3 sphericalCoord)
{
    float radius = sphericalCoord.x;
    float phi = sphericalCoord.y;
    float theta = sphericalCoord.z;
    return vec3(
        radius * sin(phi) * cos(theta),
        radius * cos(phi),
        radius * sin(phi) * sin(theta)
    );
}

vec3 getRightVector(vec3 coord)
{
    return vec3(sin(coord.z), 0.0, -cos(coord.z));
}

vec3 reflection(vec3 v, vec3 n)
{
    return -v + 2.0 * dot(v, n) * n;
}

float intensity(vec3 eye, vec3 p, vec3 n,
    float kSpec, float specWeight, float diffWeight)
{
    vec3 toLight = normalize(light1 - p);
    vec3 toEye = normalize(eye - p);
    vec3 ref = normalize(reflection(toLight, n));

    float diffuse = 0.0;
    float specular = 0.0;

    // Do diffuse lighting
    diffuse +=  max(0.0, dot(toLight, n));
    // Do specular lighting
    specular += pow(max(dot(ref, toEye), 0.0), kSpec);

    // Recompute for second light
    toLight = normalize(light2 - p);
    ref = normalize(reflection(toLight, n));
    diffuse += max(0.0, dot(toLight, ref));
    specular += pow(max(dot(ref, toEye), 0.0), kSpec);

    // Recompute for third light
    toLight = normalize(light3 - p);
    ref = normalize(reflection(toLight, n));
    diffuse += max(0.0, dot(toLight, ref));
    specular += pow(max(dot(ref, toEye), 0.0), kSpec);

    return specWeight * specular + diffWeight * diffuse;
}

vec4 shadeSphere(vec3 eye, vec3 pos, Particle p)
{
    vec3 n = normalize(pos - p.pos);
    return 0.5 * intensity(eye, pos, n, 0.0, 0.0, 1.0) * p.color;
}

vec4 shade(vec3 eye, vec3 pos, vec3 n, vec4 color)
{
    return 0.5 * intensity(eye, pos, n, 0.0, 0.0, 1.0) * color;
}

vec3 getNormal(Particle p, vec3 pos)
{
    return normalize(pos - p.pos);
}

float getFieldValue(Particle p, vec3 pos)
{
    //vec3 relPos = pos - p.pos;
    float r = length(pos - p.pos)/p.radius; // normalize r so that it's from 0 to 1
    if (r > 1.0)
        return 0.0; //out of bound
    float g = r * r * r * (r * (r * 6.0 - 15.0) + 10.0);
    return 1.0 - g; //we want g(0) = 1, g(1) = 0;
}

float getFieldValue(float R, float d)
{
    float r = d/R;
    if (r > 1.0)
        return 0.0;
    float g = r * r * r * (r * (r * 6.0 - 15.0) + 10.0);
    return 1.0 - g;
}

float frac(float x)
{
    return x - floor(x);
}

bool isOutside(Particle p, vec3 pos)
{
    return length(pos - p.pos) > p.radius;
}

float getDistance(Particle p, vec3 pos)
{
    return length(pos - p.pos);
}

void main()
{
    vec3 eye     = getCartesian(eyeCoordinate);
    vec3 focus   = vec3(0, 0, 0);
    vec3 forward = normalize(focus - eye);
    vec3 right   = normalize(getRightVector(eyeCoordinate));
    vec3 up      = normalize(cross(right, forward));

    float f = 2.0;
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;

    float ar = resolution.x / resolution.y;
    right *= ar;

    vec4 background = vec4(0.0);

    vec3 imagePos = eye + right * u + up * v + forward * f;
    vec3 dir = normalize(imagePos - eye); // perspective view

    // Temporary hard coded particles
    Particle p[PARTICLESNUM];
    p[0].pos = vec3(sin(time)+0.5,0.0,0.0);
    p[0].radius = 0.5;
    p[0].color = vec4(1.0,0.0,0.0,0.0); // Red

    p[1].pos = vec3(-sin(time)-0.5,0.0,0.0);
    p[1].radius = 1.0;
    p[1].color = vec4(0.0,1.0,0.0,0.0); // Green

    gl_FragColor = background;

    float minT = -1.0;

    float stepSize = -1.0;
    vec3 pos = imagePos;
    for (int i = 0; i < MAXSTEP; i++)
    {
        // Evaluate from each particles
        float val = 0.0;
        vec3 n;
        vec4 color;
        for (int j = 0; j < PARTICLESNUM; j++)
        {   
            vec3 normal = pos - p[j].pos;
            float d = length(normal);
            if (stepSize < 0.0 || d < stepSize)
                stepSize = d;
            if (d > p[j].radius)
                continue;
            float w = getFieldValue(p[j].radius, d);
            val += w;
            n += w*normalize(normal);
            color += w*p[j].color;
        }
        if (val > THRESHOLD)
        {
            n /= val;
            color /= val;
            gl_FragColor = shade(eye, pos, n, color);
        }

        pos += STEP;        
    }

    return;
}