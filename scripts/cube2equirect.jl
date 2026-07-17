# convert a cubemap skybox (a directory of px/nx/py/ny/pz/nz face PNGs, the
# manual_wormhole format) to the equirectangular PPM TexturedSky samples.
# Face (u,v) conventions replicate manual_wormhole src/deadsimple.rs
# Skybox::sample exactly:
#   px: u=-z v=y | nx: u=z v=y | py: u=x v=-z | ny: u=x v=z
#   pz: u=x v=y  | nz: u=-x v=y
# with pixel x=(u+1)/2*w (col from left), y=(v+1)/2*h (row from top); bilinear
# where the Rust samples nearest, clamped at face borders. Faces are read via
# `magick face.png face.ppm` into a temp dir, so ImageMagick must be on PATH.
#   julia --project --threads=auto scripts/cube2equirect.jl <face_dir> <out.ppm> [WxH]
# e.g.
#   julia --project scripts/cube2equirect.jl ~/projects/manual_wormhole/res/skyboxes/bg0 out/space1.ppm
using jupiter
const J = jupiter

length(ARGS) >= 2 || error("usage: cube2equirect.jl <face_dir> <out.ppm> [WxH]")
const face_dir = ARGS[1]
const out_path = ARGS[2]
const w, h = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], "x")) : (4096, 2048)

function face_sample(img, u, v)
    _, fw, fh = size(img)
    xf = (u + 1) / 2 * fw - 0.5 # 0-based texel-center coords
    yf = (v + 1) / 2 * fh - 0.5
    x0 = floor(Int, xf); tx = xf - x0
    y0 = floor(Int, yf); ty = yf - y0
    cx(k) = clamp(k, 0, fw - 1) + 1
    cy(k) = clamp(k, 0, fh - 1) + 1
    texel(i, j) = @view img[:, cx(i), cy(j)]
    (1 - ty) * ((1 - tx) * texel(x0, y0) + tx * texel(x0 + 1, y0)) +
    ty * ((1 - tx) * texel(x0, y0 + 1) + tx * texel(x0 + 1, y0 + 1))
end

function sample_cube(fs, d)
    a = abs.(d)
    if a[1] >= a[2] && a[1] >= a[3]
        d[1] > 0 ? face_sample(fs.px, -d[3] / a[1], d[2] / a[1]) :
                   face_sample(fs.nx, d[3] / a[1], d[2] / a[1])
    elseif a[2] >= a[3]
        d[2] > 0 ? face_sample(fs.py, d[1] / a[2], -d[3] / a[2]) :
                   face_sample(fs.ny, d[1] / a[2], d[3] / a[2])
    else
        d[3] > 0 ? face_sample(fs.pz, d[1] / a[3], d[2] / a[3]) :
                   face_sample(fs.nz, -d[1] / a[3], d[2] / a[3])
    end
end

tmp = mktempdir()
faces = (; (Symbol(f) => (run(`magick $(joinpath(face_dir, f * ".png")) $(joinpath(tmp, f * ".ppm"))`);
                          J.load_ppm(joinpath(tmp, f * ".ppm")))
            for f in ("px", "nx", "py", "ny", "pz", "nz"))...)

img = zeros(3, w, h)
Threads.@threads for j in 1:h
    for i in 1:w
        az = (i - 0.5) / w * 2pi - pi
        el = pi / 2 - (j - 0.5) / h * pi
        d = [cos(el) * cos(az), cos(el) * sin(az), sin(el)]
        img[:, i, j] = sample_cube(faces, d)
    end
end

open(out_path, "w") do io
    write(io, "P6\n$w $h\n255\n")
    for j in 1:h, i in 1:w
        write(io, UInt8.(round.(Int, clamp.(img[:, i, j], 0, 1) * 255)))
    end
end
println("wrote ", out_path, " (", w, "x", h, ")")
