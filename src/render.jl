# renderer assembly: everything here lives in ambient flat spaces; the throat
# interior is trace_geodesic's business

struct Scene{MS <: NTuple{2, Mouth}, S} # one mouth per side for now; each side's ambient space is its own universe
    throat::Throat
    mouths::MS # parametric so per-passage mouth dispatch stays static
    sky::S # sky(side, dir) -> rgb in [0,1]^3
end

"""
compute budget per camera ray. Wormholes admit limit cycles — rays that loop
through the throat(s) forever instead of escaping — so whoever renders decides
how hard to chase near-limit rays: step_size and max_steps bound one throat
passage, max_passages bounds how many times a ray may re-enter a mouth.
tolerance selects the integrator: 0 (the default) is fixed-step RK4 at
step_size — the cheaper choice at image accuracy on fine-chart meshes, where
the per-step chart-displacement cap leaves the DP5 pair no room to stretch;
tolerance > 0 is the error-controlled Dormand–Prince tracer (step_size seeds
it, max_steps counts attempts), which buys explicit accuracy and wins on
coarse-chart scenes or at tight tolerances.
"""
struct RayBudget
    step_size::Float64
    max_steps::Int
    max_passages::Int
    tolerance::Float64
end

RayBudget(step_size, max_steps, max_passages) = RayBudget(step_size, max_steps, max_passages, 0.0)

"""
follow an ambient ray through however many throat passages the budget allows;
returns the escaping AmbientRay, or nothing if the budget ran out first
"""
function trace_ray(env, scene, budget, ray::AmbientRay)
    for _ in 1:budget.max_passages
        mouth = scene.mouths[side(half_throat(ray))]
        v = enter_mouth(env, mouth, ray.pos, ray.vel)
        v === nothing && return ray # nothing left to hit but sky
        res = budget.tolerance > 0 ?
              trace_geodesic(env, v, budget.step_size, budget.max_steps, budget.tolerance) :
              trace_geodesic(env, v, budget.step_size, budget.max_steps)
        res isa SituatedPhase && return nothing # out of steps inside the throat
        ray = res
    end
    nothing # out of passages; presumably orbiting near a limit cycle
end

"""
pinhole camera: a point plus a frame (columns: right, up, forward). The frame
is any basis, NOT necessarily orthonormal — orthonormality is merely
look_at_camera's choice of initial condition. Transporting a camera through
wormholes is connection transport, not metric transport (with non-isometric
placements there is no global metric), so a camera returning from a loop may
come back rescaled or sheared; its images will honestly show that. A uniform
rescale of the frame is invisible in a single image: ray directions are homogeneous in it.
However, the scale of the frame may be observed in a camera transport process,
if the velocity is proportional to the frame.
"""
struct Camera
    pos::SVector{3, Float64}
    frame::SMatrix{3, 3, Float64, 9} # any basis, see above — concrete type, not orthonormality, is what's fixed here
    tan_half_fov::Float64
end

function look_at_camera(pos, target, up, fov)
    forward = generic_normalize(target - pos)
    right = generic_normalize(cross(forward, up))
    Camera(pos, [right cross(right, forward) forward], tan(fov / 2))
end

function camera_ray(camera::Camera, half_throat::HalfThroat, x, y) # x, y in [-1, 1], x aspect-scaled
    dir = generic_normalize(camera.frame * SVector(x * camera.tan_half_fov, y * camera.tan_half_fov, 1.0))
    AmbientRay(half_throat, camera.pos, dir)
end

function pixel_ray(camera::Camera, half_throat::HalfThroat, i, j, width, height) # pixel centers, row 1 at the top
    x = (2 * (i - 0.5) / width - 1) * (width / height)
    y = 1 - 2 * (j - 0.5) / height
    camera_ray(camera, half_throat, x, y)
end

function checker_sky(side, dir)
    az = atan(dir[2], dir[1])
    el = asin(clamp(dir[3], -1.0, 1.0))
    parity = mod(floor(Int, (az + pi) / (pi / 6)) + floor(Int, (el + pi / 2) / (pi / 6)), 2)
    hi = side == 1 ? [0.35, 0.55, 0.95] : [0.95, 0.5, 0.2]
    lo = side == 1 ? [0.08, 0.12, 0.3] : [0.35, 0.12, 0.05]
    parity == 0 ? hi : lo
end

unresolved_color() = [0.0, 0.0, 0.0]

struct RayMap # per-pixel trace outcome; trace once, shade under as many skies as you like
    side::Matrix{Int} # exit side; 0 = unresolved within budget
    pos::Array{Float64, 3} # 3 x width x height, ambient exit positions
    vel::Array{Float64, 3} # 3 x width x height, ambient exit velocities
end

"""
the expensive pass: trace every camera ray through the scene and record where
it ends up, without committing to any particular sky
"""
function render_raymap(env, scene, budget, camera::Camera, cam_side::Int, width::Int, height::Int)
    raymap = RayMap(zeros(Int, width, height), zeros(3, width, height), zeros(3, width, height))
    cam_space = HalfThroat(scene.throat, cam_side)
    Threads.@threads for j in 1:height
        for i in 1:width
            out = trace_ray(env, scene, budget, pixel_ray(camera, cam_space, i, j, width, height))
            out === nothing && continue
            raymap.side[i, j] = side(half_throat(out))
            raymap.pos[:, i, j] = out.pos
            raymap.vel[:, i, j] = out.vel
        end
    end
    raymap
end

"""
the cheap pass: rgb image (3 x width x height, values in [0,1]) from a raymap
"""
function shade(raymap::RayMap, sky)
    w, h = size(raymap.side)
    img = zeros(3, w, h)
    for j in 1:h, i in 1:w
        img[:, i, j] = raymap.side[i, j] == 0 ? unresolved_color() :
                       sky(raymap.side[i, j], generic_normalize(raymap.vel[:, i, j]))
    end
    img
end

function render(env, scene, budget, camera::Camera, cam_side::Int, width::Int, height::Int)
    shade(render_raymap(env, scene, budget, camera, cam_side, width, height), scene.sky)
end

"""
straight-ray reference pass: headlight-shade the camera-side mouth
tessellation (grayscale, silhouette-darkening |n·view|), sky everywhere else —
the undeflected outline to hold lensed renders against
"""
function render_flat(env, scene, camera::Camera, cam_side::Int, width::Int, height::Int)
    img = zeros(3, width, height)
    cam_space = HalfThroat(scene.throat, cam_side)
    mouth = scene.mouths[cam_side]
    Threads.@threads for j in 1:height
        for i in 1:width
            ray = pixel_ray(camera, cam_space, i, j, width, height)
            hit = nearest_mouth_hit(mouth, ray.pos, ray.vel)
            if hit === nothing
                img[:, i, j] = scene.sky(cam_side, ray.vel)
            else
                _, tri = hit
                n = generic_normalize(cross(tri.corners[2] - tri.corners[1], tri.corners[3] - tri.corners[1]))
                img[:, i, j] .= 0.2 + 0.8 * abs(n' * ray.vel)
            end
        end
    end
    img
end

save_raymap(path, raymap::RayMap) = open(io -> serialize(io, raymap), path, "w")
load_raymap(path) = open(deserialize, path)

"""
equirectangular texture lookup: azimuth wraps around the image width, elevation
runs top row = zenith; bilinear filtering, azimuth-wrapping and pole-clamping.
img is 3 x width x height in [0,1], the layout save_ppm/load_ppm use.
"""
function sample_equirect(img, dir)
    _, w, h = size(img)
    az = atan(dir[2], dir[1])
    el = asin(clamp(dir[3], -1.0, 1.0))
    x = (az + pi) / (2pi) * w + 0.5 # continuous pixel coords, texel centers at integers
    y = (1 - (el + pi / 2) / pi) * h + 0.5
    x0 = floor(Int, x); tx = x - x0
    y0 = floor(Int, y); ty = y - y0
    xi(k) = mod(k - 1, w) + 1
    yi(k) = clamp(k, 1, h)
    texel(i, j) = SVector(img[1, i, j], img[2, i, j], img[3, i, j])
    (1 - ty) * ((1 - tx) * texel(xi(x0), yi(y0)) + tx * texel(xi(x0 + 1), yi(y0))) +
    ty * ((1 - tx) * texel(xi(x0), yi(y0 + 1)) + tx * texel(xi(x0 + 1), yi(y0 + 1)))
end

struct TexturedSky # a Scene-ready sky: one equirectangular texture per side
    images::NTuple{2, Array{Float64, 3}}
end

TexturedSky(img1, img2) = TexturedSky((img1, img2))

(sky::TexturedSky)(side, dir) = sample_equirect(sky.images[side], dir)

function save_ppm(path, img)
    _, w, h = size(img)
    open(path, "w") do io
        println(io, "P3\n", w, " ", h, "\n255")
        for j in 1:h, i in 1:w
            println(io, join(round.(Int, clamp.(img[:, i, j], 0, 1) * 255), " "))
        end
    end
end

function ppm_token(io) # next whitespace-delimited header token; '#' comments run to end of line
    tok = IOBuffer()
    while !eof(io)
        ch = read(io, Char)
        if ch == '#'
            readline(io)
        elseif isspace(ch)
            position(tok) > 0 && return String(take!(tok))
        else
            write(tok, ch)
        end
    end
    String(take!(tok))
end

"""
reads a P3 (ASCII) or P6 (binary, 8-bit) PPM into the 3 x width x height [0,1]
layout save_ppm writes. The route for real images: `magick x.png x.ppm` emits P6.
"""
function load_ppm(path)
    open(path) do io
        magic = ppm_token(io)
        w = parse(Int, ppm_token(io))
        h = parse(Int, ppm_token(io))
        maxval = parse(Int, ppm_token(io))
        img = zeros(3, w, h)
        if magic == "P6"
            maxval <= 255 || error("16-bit P6 unsupported (maxval $maxval)")
            data = read(io, 3 * w * h) # binary payload starts right after the single whitespace ppm_token consumed
            length(data) == 3 * w * h || error("truncated P6 payload in $path")
            for j in 1:h, i in 1:w, c in 1:3
                img[c, i, j] = data[3 * ((j - 1) * w + i - 1) + c] / maxval
            end
        elseif magic == "P3"
            for j in 1:h, i in 1:w, c in 1:3
                img[c, i, j] = parse(Int, ppm_token(io)) / maxval
            end
        else
            error("unsupported PPM magic '$magic' in $path")
        end
        img
    end
end
