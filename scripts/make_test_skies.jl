# generates the two tracked test skies in res/textures/: equirectangular
# graticule patterns (white lines every 30 degrees) over a cool palette for
# side 1 and a warm one for side 2 — hue tracks azimuth, brightness elevation,
# so orientation flips and seam errors are visible at a glance. Rerun to
# regenerate; real panoramas come in via `magick photo.png sky.ppm` instead.
using jupiter
const J = jupiter

function test_sky(w, h, warm::Bool)
    img = zeros(3, w, h)
    for j in 1:h, i in 1:w
        az = (i - 0.5) / w * 360 - 180
        el = 90 - (j - 0.5) / h * 180
        t = (el + 90) / 180
        s = 0.5 + 0.5 * sind(az)
        cool = [0.08 + 0.15 * (1 - t), 0.12 + 0.5 * t * s, 0.3 + 0.6 * t]
        heat = [0.3 + 0.6 * t, 0.12 + 0.5 * t * s, 0.08 + 0.15 * (1 - t)]
        on_grid = any(abs(mod(x + 15, 30) - 15) < 0.6 for x in (az, el))
        img[:, i, j] = on_grid ? [1.0, 1.0, 1.0] : (warm ? heat : cool)
    end
    img
end

function save_p6(path, img)
    _, w, h = size(img)
    open(path, "w") do io
        write(io, "P6\n$w $h\n255\n")
        for j in 1:h, i in 1:w
            write(io, UInt8.(round.(Int, clamp.(img[:, i, j], 0, 1) * 255)))
        end
    end
end

dir = joinpath(@__DIR__, "..", "res", "textures")
mkpath(dir)
# 512x256 keeps each file under jj's 1 MiB new-file guard; ample for a test pattern
save_p6(joinpath(dir, "sky1.ppm"), test_sky(512, 256, false))
save_p6(joinpath(dir, "sky2.ppm"), test_sky(512, 256, true))
println("wrote ", joinpath(dir, "sky1.ppm"), " and sky2.ppm (512x256 P6)")
