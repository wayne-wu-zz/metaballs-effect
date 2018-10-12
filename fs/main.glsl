#extension GL_OES_standard_derivatives : enable
#ifdef GL_ES
precision highp float;
#endif

uniform float time;
uniform vec2 resolution;
uniform vec3 eyeCoordinate;

uniform int texWidth;
uniform int texHeight;
uniform int particleNum;
uniform sampler2D positionTexture;

uniform vec3 light1;
uniform vec3 light2;
uniform vec3 light3;

bool sphereIntersections(vec3 c, float r, vec3 o, vec3 d, inout float t1, inout float t2)
{
    // Ray/Sphere Intersection Test
    // float A = 1;
    vec3 v = o - c;
    float B = 2.0 * dot(d, v);
    float C = dot(v, v) - r*r;

    float det = B*B - 4.0*C;

    if (det <= 0.0)
        return false; // No solution

    t1 = -B - sqrt(det);
    t2 = -B + sqrt(det);

    if (t1 < 0.0 || t2 < 0.0)
        return false; // Sphere is behind the ray

    return true; 
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

float frac(float a, float b)
{
    return a - b * floor(a/b);
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

vec4 shade(vec3 eye, vec3 pos, vec3 n, vec4 color)
{
    return 0.5 * intensity(eye, pos, n, 0.0, 0.0, 1.0) * color;
}

/* Get the texture coordinate for the given particle index */
vec2 getTexCoord(int idx)
{
    float i = float(idx);
    float w = float(texWidth);
    float h = float(texHeight);
    float y = floor(i/h) + 0.5;
    float x = frac(i,h)/w + 0.5;
    return vec2(x/w, y/h);
}

/* Signed Distance Field for Sphere */
float sdSphere( vec3 pos, vec3 center, float r)
{
    return length(pos - center) - r;
}

/* Get the field value based on the radius */
float getFieldValue(float R, float d)
{
    float r = d/R;
    if (r > 1.0)
        return 0.0;
    float g = r * r * r * (r * (r * 6.0 - 15.0) + 10.0);
    return 1.0 - g;
}

float sdMetaballs( vec3 pos )
{
    return 0.0;
}

void main()
{
    vec3 eye     = getCartesian(eyeCoordinate);
    vec3 focus   = vec3(0, 0, 0);
    vec3 forward = normalize(focus - eye);
    vec3 right   = normalize(getRightVector(eyeCoordinate));
    vec3 up      = normalize(cross(right, forward));

    float f = 1.0;
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;

    float ar = resolution.x / resolution.y;
    right *= ar;

    vec4 background = vec4(0.0);

    vec3 imagePos = eye + right * u + up * v + forward * f;

    vec3 dir = normalize(imagePos - eye);

    float minstep = 10000.0;
    int minparticle = -1;

    vec3 pos = imagePos;
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            vec4 v = texture2D(positionTexture, getTexCoord(j));

            float sdf = sdSphere(pos, v.xyz, v.w);
            if (sdf < minstep)
            {
                minstep = sdf;
                minparticle = j;
            }
        }

        if (minstep < 0.001)
        {
            vec4 v = texture2D(positionTexture, getTexCoord(minparticle));
            gl_FragColor = shade(eye, pos, normalize(pos-v.xyz), vec4(1.0,0.0,0.0,0.0));;
            return;
        }

        pos += minstep*dir;
    }


    gl_FragColor = background;

    return;
}