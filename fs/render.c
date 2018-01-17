#define PARTICLESNUM 2
#define THRESHOLD 0.01
#define DT 0.1
#define STEP 0.01
#define ENTRY 0
#define EXIT 1

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

/*
// Iterative Merge Sort
void mergeSort(Intersection[] arr, int n)
{
   int curr_size;  // For current size of subarrays to be merged
                   // curr_size varies from 1 to n/2
   int left_start; // For picking starting index of left subarray
                   // to be merged
 
   // Merge subarrays in bottom up manner.  First merge subarrays of
   // size 1 to create sorted subarrays of size 2, then merge subarrays
   // of size 2 to create sorted subarrays of size 4, and so on.
   for (curr_size=1; curr_size<=n-1; curr_size = 2*curr_size)
   {
       // Pick starting point of different subarrays of current size
       for (left_start=0; left_start<n-1; left_start += 2*curr_size)
       {
           // Find ending point of left subarray. mid+1 is starting 
           // point of right
           int mid = left_start + curr_size - 1;
 
           int right_end = min(left_start + 2*curr_size - 1, n-1);
 
           // Merge Subarrays arr[left_start...mid] & arr[mid+1...right_end]
           merge(arr, left_start, mid, right_end);
       }
   }
}

*/

/*
void merge(Intersection arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 =  r - m;
 
    // create temp arrays 
    Intersection L[n1], R[n2];
 
    // Copy data to temp arrays L[] and R[]
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1+ j];
 
    // Merge the temp arrays back into arr[l..r]
    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2)
    {
        if (L[i].t <= R[j].t)
        {
            arr[k] = L[i];
            i++;
        }
        else
        {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    // Copy the remaining elements of L[], if there are any
    while (i < n1)
    {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    // Copy the remaining elements of R[], if there are any
    while (j < n2)
    {
        arr[k] = R[j];
        j++;
        k++;
    }
}
*/

vec3 getNormal(Particle p, vec3 pos)
{
    return normalize(pos - p.pos);
}

float getFieldValue(Particle p, vec3 pos)
{
    vec3 relPos = pos - p.pos;
    float r = p.radius;
    float g = r * r * r * (r * (r * 6.0 - 15.0) + 10.0);
    return g;
}

float frac(float x)
{
    return x - floor(x);
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

    Particle p[PARTICLESNUM];
    p[0].pos = vec3(sin(time)+1.0,0.0,0.0);
    p[0].radius = 1.5;
    p[0].color = vec4(1.0,0.0,0.0,0.0); // Red

    p[1].pos = vec3(-sin(time)-1.0,0.0,0.0);
    p[1].radius = 1.5;
    p[1].color = vec4(0.0,1.0,0.0,0.0); // Green

    gl_FragColor = background;

    // Find all the intersections with the particles
    const int maxNumIntersections = PARTICLESNUM*2;
    float intersections[maxNumIntersections];
    Intersection arr[maxNumIntersections];
    float minT = -1.0;
    Particle minP;
    for (int i = 0; i < 2; i++)
    {
        float t1, t2;
        bool hitSphere = sphereIntersections(p[i], imagePos, dir, t1, t2);
        if (hitSphere)
        {
            if (minT < 0.0 || min(t1, t2) < minT)
            {
                minT = min(t1, t2);
                minP = p[i];
            }
            intersections[i] = t1;
            intersections[i+1] = t2;
            //arr[i] = Intersection(t1, i);
            //arr[i+1] = Intersection(t2, i)
        }
        else
        {
            // Set them to -1.0 indicating that there's no intersection
            intersections[i] = -1.0;
            intersections[i+1] = -1.0;
            //arr[i] = Intersection(-1.0, i);
            //arr[i+1] = Intersection(-1.0, i);
        }
    }

    if (minT < 0.0)
    {
        return; // no intersection
    }

    /*
    mergeSort(arr, maxNumIntersections);

    Particle active[PARTICLESNUM];
    int active_size = 0;
    for (int i = 0; i < maxNumIntersections; i++)
    {
        Intersection tmp = arr[i];
        if (tmp.type == ENTRY)
            p[tmp.pIdx].active = true;
        else if (tmp.type == EXIT)
            p[tmp.pIdx].active = false;

        float pos = imagePos + dir*tmp.t;
        while(true)
        {
            pos += STEP;
        }
        for (int j = 0; j < PARTICLESNUM; j++)
        {
            if (p[j].active)
            {
                
            }
        }

    }
    */

    // float totalWeight = 0.0;
    // vec3 normal;
    // vec4 color;
    // vec3 pos;
    // for (int i = 0; i < maxNumIntersections; i++)
    // {
    //     float t = intersections[i];
    //     if (t < 0.0)
    //         continue; // no intersection

    //     float dt = t - minT;
    //     if (dt < DT)
    //     {
    //         float w = 1.0 - dt/DT;
    //         totalWeight += w;
    //         vec3 hitPos = imagePos + dir*t;
    //         pos += w*hitPos;
    //         normal += w*getNormal(p[i/2], hitPos);
    //         color += w*p[i/2].color;
    //         //gl_FragColor += w*shadeSphere(eye, imagePos + dir*t, p[i/2]);
    //     } 
    // }
    // pos /= totalWeight;
    // normal /= totalWeight;
    // color /= totalWeight;
    // gl_FragColor = shade(eye, pos, normal, color);
    

    gl_FragColor = shadeSphere(eye, imagePos + dir*minT, minP);
}