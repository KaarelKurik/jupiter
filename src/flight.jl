# piecewise-geodesic camera flight: a FlyingCamera coasts along geodesics
# (parallel-transporting its frame), crossing mouths in either direction with
# a continuous geodesic parameter; steering happens only at explicit
# frame-relative maneuvers between coasts. Choreography code — runs once per
# leg per keyframe, so clarity over allocation discipline here.

"""
a camera in flight: state is an AmbientRay (frame in ambient coordinates) or
a SituatedPhase (frame in chart coordinates); the state's vel is the motion
"""
struct FlyingCamera{S, F}
    state::S
    frame::F
end

with_vel(s::AmbientRay, vel) = AmbientRay(s.half_throat, s.pos, vel)
with_vel(s::SituatedPhase, vel) = SituatedPhase(s.chart, s.pos, vel)

"""
frame-relative re-orientation: the new frame is E * R, i.e. R's columns are
the new axes expressed in the camera's own current basis
"""
maneuver(cam::FlyingCamera, R) = FlyingCamera(cam.state, cam.frame * R)

"""
frame-relative re-aim of the motion: vel becomes frame * u (u in frame
coordinates), e.g. steer(cam, [0, 0, s]) coasts forward at speed s
"""
steer(cam::FlyingCamera, u) = FlyingCamera(with_vel(cam.state, cam.frame * u), cam.frame)

"""
advance the camera by geodesic parameter Δτ (in units of its vel), stepping at
h inside the throat, crossing mouths in either direction as needed; the
parameter is continuous across crossings (entry pullback and collar export
both preserve vel). max_crossings guards grazing pathologies.
"""
function coast(env, scene, cam::FlyingCamera, Δτ, h; max_crossings=8)
    state, E = cam.state, cam.frame
    τ = float(Δτ)
    ε = 1e-12 # don't chase float-remainder micro-steps (0.5 is not 25 exact steps of 0.02)
    for _ in 1:max_crossings
        τ > ε || break
        if state isa AmbientRay
            entry = enter_transport(env, scene.mouths[side(half_throat(state))], state.pos, state.vel, E)
            if entry === nothing || entry[3] > τ
                state = AmbientRay(half_throat(state), state.pos + τ * state.vel, state.vel)
                τ = 0.0
            else
                state, E = entry[1], entry[2]
                τ -= entry[3]
            end
        else
            while τ > ε
                state, E = settle_transport(env, transport_step(env, state, E, min(h, τ))...)
                τ -= min(h, τ)
                if exits_mouth(state)
                    state, E = to_ambient(env, state, E)
                    break
                end
            end
        end
    end
    FlyingCamera(state, E)
end

"""
render the flying camera's view: ambient states through the pinhole Camera,
in-chart states through SituatedCamera emission — the continuity of the two
descriptions across a crossing is certified in the tests
"""
function keyframe_raymap(env, scene, budget, cam::FlyingCamera, tan_half_fov, width::Int, height::Int)
    if cam.state isa AmbientRay
        camera = Camera(cam.state.pos, cam.frame, tan_half_fov)
        render_raymap(env, scene, budget, camera, side(half_throat(cam.state)), width, height)
    else
        render_raymap(env, scene, budget, SituatedCamera(cam.state.chart, cam.state.pos, cam.frame, tan_half_fov), width, height)
    end
end
