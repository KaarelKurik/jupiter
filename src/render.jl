# renderer assembly: everything here lives in ambient flat spaces; the throat
# interior is trace_geodesic's business

struct Scene # one mouth per side for now; each side's ambient space is its own universe
    throat::Throat
    mouths::NTuple{2, Mouth}
    sky # sky(side, dir) -> rgb in [0,1]^3
end

"""
compute budget per camera ray. Wormholes admit limit cycles — rays that loop
through the throat(s) forever instead of escaping — so whoever renders decides
how hard to chase near-limit rays: step_size and max_steps bound one throat
passage, max_passages bounds how many times a ray may re-enter a mouth.
"""
struct RayBudget
    step_size
    max_steps
    max_passages
end

"""
follow an ambient ray through however many throat passages the budget allows;
returns the escaping AmbientRay, or nothing if the budget ran out first
"""
function trace_ray(env, scene, budget, ray::AmbientRay)
    for _ in 1:budget.max_passages
        mouth = scene.mouths[side(half_throat(ray))]
        v = enter_mouth(env, mouth, ray.pos, ray.vel)
        v === nothing && return ray # nothing left to hit but sky
        res = trace_geodesic(env, v, budget.step_size, budget.max_steps)
        res isa SituatedPhase && return nothing # out of steps inside the throat
        ray = res
    end
    nothing # out of passages; presumably orbiting near a limit cycle
end

struct Camera
    pos
    frame # columns: right, up, forward
    tan_half_fov
end

function look_at_camera(pos, target, up, fov)
    forward = generic_normalize(target - pos)
    right = generic_normalize(cross(forward, up))
    Camera(pos, [right cross(right, forward) forward], tan(fov / 2))
end

function camera_ray(camera::Camera, half_throat::HalfThroat, x, y) # x, y in [-1, 1], x aspect-scaled
    dir = generic_normalize(camera.frame * [x * camera.tan_half_fov, y * camera.tan_half_fov, 1.0])
    AmbientRay(half_throat, camera.pos, dir)
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

"""
render the scene from a camera sitting in cam_side's ambient space;
returns a 3 x width x height array of rgb values in [0,1]
"""
function render(env, scene, budget, camera::Camera, cam_side::Int, width::Int, height::Int)
    img = zeros(3, width, height)
    cam_space = HalfThroat(scene.throat, cam_side)
    aspect = width / height
    for j in 1:height, i in 1:width
        x = (2 * (i - 0.5) / width - 1) * aspect
        y = 1 - 2 * (j - 0.5) / height
        out = trace_ray(env, scene, budget, camera_ray(camera, cam_space, x, y))
        img[:, i, j] = out === nothing ? unresolved_color() :
                       scene.sky(side(half_throat(out)), generic_normalize(out.vel))
    end
    img
end

function save_ppm(path, img)
    _, w, h = size(img)
    open(path, "w") do io
        println(io, "P3\n", w, " ", h, "\n255")
        for j in 1:h, i in 1:w
            println(io, join(round.(Int, clamp.(img[:, i, j], 0, 1) * 255), " "))
        end
    end
end
