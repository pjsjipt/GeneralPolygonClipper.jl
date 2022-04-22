module GeneralPolygonClipper


@enum GPCOperation GPC_DIFF=0 GPC_INT GPC_XOR GPC_UNION


export Vertex, area, centroid
export gpc_polygon_clip, gpc_polygon_to_tristrip, gpc_tristrip_clip
export GPCPolygon, GPCTriStrip
export union_strip, intersect_strip, diff_strip, xor_strip
export gpc2tristrip, pintri
export GPCOperation, GPC_DIFF, GPC_INT, GPC_XOR, GPC_UNION


struct Vertex
    x::Float64
    y::Float64
end

Vertex(v) = Vertex(v[1], v[2])


"""
`gpc_vertex_list(nv, vertex)`

An `struct` to be passed to GPC C code library that stores a vector of
vertices.
"""
struct gpc_vertex_list
    "Number of vertices"
    num_vertices::Cint
    "Pointer to the vector of vertices"
    vertex::Ptr{Vertex}
end


"""
`gpc_polygon(nc, hole, contour)`

An `struct` to be passed to GPC C code library that stores a polygon.
In the GPC library, a polygon is a sequence of contours. 
Each contour is a sequence of vertices and it can correspond to a hole or not.

"""
struct gpc_polygon
    "Number of contours of the polygon"
    num_contours::Cint
    "Vector specifying whether a contour is a hole or not"
    hole::Ptr{Cint}
    "Vector of [`gpc_vertex_list`](@ref)"
    contour::Ptr{gpc_vertex_list}
end

"""
`GPCPolygon(holes, contours)`

A GPC Polygon where each polygon is composed of 1 or more contours.
Each contour is a vector of vertices and can be a hole or not.
"""
struct GPCPolygon
    "Is the i-th contour a hole?"
    holes::Vector{Bool}
    "Vector of vertices of the i-th contour"
    contours::Vector{Vector{Vertex}}
end

function GPCPolygon(x::AbstractVector, y::AbstractVector)
    nv = length(x)
    @assert nv == length(y)

    v = [Vertex(x[i], y[i]) for i in 1:nv]

    return GPCPolygon([false], [v])
end



"""
`area(p::GPCPolygon)`

Calculates the surface area of a [`GPCPolygon`](@ref) polygon.

The area is defined as the sum of the areas of external contours minus
the area of the holes. 

There is an issue here as to the signal of the area. 
Depending on how the vertices on a contour are oriented, 
the area can be negative. To keep things consistent, the area 
of external contours are positive and the area of holes is negative.

"""
function area(p::GPCPolygon)

    A = 0.0

    for (i, c) in enumerate(p.contours)
        if p.holes[i]
            A -= abs(area(c))
        else
            A += abs(area(c))
        end
    end
    return A
    
end

"""
`area(c::Vector{Vertex})`

Compute the area of a closed contour of vertices.
"""
function area(c::Vector{Vertex})
    A = c[end].x * c[1].y - c[1].x*c[end].y
    nv = length(c)
    for i in 2:nv
        A = A + c[i-1].x * c[i].y - c[i].x * c[i-1].y
    end

    return A/2
    
end

"""
`triarea(u, v, w)`

Calculates the area of the triangle defined by vertices `u`, `v`, `w`.
"""
triarea(u, v, w) = abs(u.x*v.y - v.x*u.y + v.x*w.y - w.x*v.y + w.x*u.y - u.x*w.y)/2

"""
`centroid(p::GPCPolygon)`

Calculates the centroid of a [`GPCPolygon`](@ref).

The function decomposes the polygon into triangle strips and uses
that to compute the centroid.
"""
function centroid(p::GPCPolygon)
    Cx = 0.0
    Cy = 0.0
    A = 0.0
    strip = gpc2tristrip(p)
    for s in strip
        for i in 1:length(s)
            v = s[i]
            Ai = triarea(v...)
            Cx += Ai * (v[1].x + v[2].x + v[3].x)
            Cy += Ai * (v[1].y + v[2].y + v[3].y)
            A += Ai
        end
    end
    
    return Vertex(Cx/3A, Cy/3A)
end

"""
`gpc_tristrip(ns, strip)`

An `struct` to be passed to GPC C code library that stores a triangle strip.
In the GPC library, a triangle strip is a sequence of `Vertex` vectors. 
Each contour is a sequence of vertices and it can correspond to a hole or not.

"""
struct gpc_tristrip
    num_strips::Cint
    strip::Ptr{gpc_vertex_list}
end

"""
`GPCStrip(strip)`

Stores triangle strips. Triangle the triangle strips follows the conventio
used by the General Polygon Clipper C library. 
"""
struct GPCTriStrip
    verts::Vector{Vertex}
end

import Base.length
length(strip::GPCTriStrip) = length(strip.verts) - 2


    
"Number of triangles in a triangle strip"
numtriangles(strip::AbstractVector{Vertex}) = length(strip.verts) - 2


const idx_odd = (1,2,3)
const idx_even = (3,2,1)


"""
`trianglevertices(ivert)`

Return the index of `ivert`-th triangle in triangle strip.
"""
trianglevertices(ivert) = 
    if ivert % 2 != 0
        return (ivert-1) .+ idx_odd
    else
        return (ivert-1) .+ idx_even
    end

import Base.getindex

function getindex(stri::GPCTriStrip, i::Integer)
    idx = trianglevertices(i)
    return stri.verts[idx[1]], stri.verts[idx[2]], stri.verts[idx[3]]
end


    

"""
`gpc_polygon_clip(op, p1, p2)`

Use GPC C library to perform boolean operations with [`GPCPolygon`](@ref).

## Arguments

 * `op` a [`GPCOperation`](@ref) specifying which operation to perform
 * `p1` first [`GPCPolygon`](@ref) 
 * `p2` second [`GPCPolygon`](@ref) 

## Returns

Returns a [`GPCPolygon`](@ref) object.
"""
function gpc_polygon_clip(op::GPCOperation, p1::GPCPolygon, p2::GPCPolygon)

    # Create corresponding C structures:
    nc1 = Cint(length(p1.holes))
    h1 = [Cint(h) for h in p1.holes]
    c1 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p1.contours]
    gpc1 = gpc_polygon(nc1, pointer(h1), pointer(c1))
                       
    nc2 = Cint(length(p2.holes))
    h2 = [Cint(h) for h in p2.holes]
    c2 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p2.contours]
    gpc2 = gpc_polygon(nc2, pointer(h2), pointer(c2))

    gpc_out = Ref(gpc_polygon(0, 0, 0))
    ccall( (:gpc_polygon_clip, "/usr/local/lib/libgpc.so"), Cvoid,
           (Cint, Ptr{gpc_polygon}, Ptr{gpc_polygon}, Ptr{gpc_polygon}),
           Cint(op), Ref(gpc1), Ref(gpc2), gpc_out)
    h = [Bool(unsafe_load(gpc_out[].hole, i)) for i in 1:gpc_out[].num_contours]
    nc = length(h)
    contours = Vector{Vertex}[]
    for i in 1:nc
        gvl = unsafe_load(gpc_out[].contour, i)
        verts = [unsafe_load(gvl.vertex, k) for k in 1:gvl.num_vertices]
        push!(contours, verts)
    end
    
    
    # Deallocate the memory allocated by the C function:
    ccall(  (:gpc_free_polygon, "/usr/local/lib/libgpc.so") , Cvoid,
            (Ptr{gpc_polygon},), gpc_out)
    
    return GPCPolygon(h, contours)
    
end

import Base.union
union(p1::GPCPolygon, p2::GPCPolygon) =  gpc_polygon_clip(GPC_UNION, p1, p2)
import Base.intersect
intersect(p1::GPCPolygon, p2::GPCPolygon) =  gpc_polygon_clip(GPC_INT, p1, p2)

import Base.(-)
(-)(p1::GPCPolygon, p2::GPCPolygon) =  gpc_polygon_clip(GPC_DIFF, p1, p2)

import Base.xor
xor(p1::GPCPolygon, p2::GPCPolygon) =  gpc_polygon_clip(GPC_XOR, p1, p2)


    
"""
`gpc_polygon_to_tristrip(p::GPCPolygon)`

`gpc2tristrip(p::GPCPolygon)`

Converts a [`GPCPolygon`](@ref) to a triangle strip.

In GPC, a triangle strip is a vector of vectors of vertices that specify a triangle strip. Actually, each vector of vertices correspond to a contiguous triangle strip and a polygon can be decomposed into several triangle strip. 

Thus this function will return a vector of vector of vertices.

In each vector of vertices, the 3 first vertices make up a triangle. Each additional
vertex adds a new triangle to the strip. See function [`trianglevertices`](@ref) to see how each triangle is made up

"""
function gpc_polygon_to_tristrip(p::GPCPolygon)
    nc1 = Cint(length(p.holes))
    h1 = [Cint(h) for h in p.holes]
    c1 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p.contours]
    gpc1 = gpc_polygon(nc1, pointer(h1), pointer(c1))

    tri_out = Ref(gpc_tristrip(0, 0))

    ccall( (:gpc_polygon_to_tristrip, "/usr/local/lib/libgpc.so"), Cvoid,
           (Ptr{gpc_polygon}, Ptr{gpc_tristrip}), Ref(gpc1), tri_out)

    ns = tri_out[].num_strips
    strip = Vector{Vertex}[]
    for i in 1:ns
        gvl = unsafe_load(tri_out[].strip, i)
        verts = [unsafe_load(gvl.vertex, k) for k in 1:gvl.num_vertices]
        push!(strip, verts)
    end

    ccall(  (:gpc_free_tristrip, "/usr/local/lib/libgpc.so") , Cvoid,
            (Ptr{gpc_tristrip},), tri_out)

    return [GPCTriStrip(v) for v in strip]
    
end


"""
`gpc_polygon_to_tristrip(p::GPCPolygon)`

`gpc2tristrip(p::GPCPolygon)`

Converts a [`GPCPolygon`](@ref) to a triangle strip.

In GPC, a triangle strip is a vector of vectors of vertices that specify a triangle strip. Actually, each vector of vertices correspond to a contiguous triangle strip and a polygon can be decomposed into several triangle strip. 

Thus this function will return a vector of vector of vertices.

In each vector of vertices, the 3 first vertices make up a triangle. Each additional
vertex adds a new triangle to the strip. See function [`trianglevertices`](@ref) to see how each triangle is made up

"""
gpc2tristrip(p::GPCPolygon) = gpc_polygon_to_tristrip(p)



"""
`gpc_tristrip_clip(op, p1, p2)`

Use GPC C library to perform boolean operations with [`GPCPolygon`](@ref).
In contrast to function [`gpc_polygon_clip`](@ref) which returns a polygon
([`GPCPolygon`](@ref)), this function returns directlty triangle strips.

## Arguments

 * `op` a [`GPCOperation`](@ref) specifying which operation to perform
 * `p1` first [`GPCPolygon`](@ref) 
 * `p2` second [`GPCPolygon`](@ref) 

## Returns

Returns a [`GPCPolygon`](@ref) object.
"""
function gpc_tristrip_clip(op::GPCOperation, p1::GPCPolygon, p2::GPCPolygon)

    # Create corresponding C structures:
    nc1 = Cint(length(p1.holes))
    h1 = [Cint(h) for h in p1.holes]
    c1 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p1.contours]
    gpc1 = gpc_polygon(nc1, pointer(h1), pointer(c1))
                       
    nc2 = Cint(length(p2.holes))
    h2 = [Cint(h) for h in p2.holes]
    c2 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p2.contours]
    gpc2 = gpc_polygon(nc2, pointer(h2), pointer(c2))

    tri_out = Ref(gpc_tristrip(0, 0))
    ccall( (:gpc_tristrip_clip, "/usr/local/lib/libgpc.so"), Cvoid,
           (Cint, Ptr{gpc_polygon}, Ptr{gpc_polygon}, Ptr{gpc_polygon}),
           Cint(op), Ref(gpc1), Ref(gpc2), tri_out)

    ns = tri_out[].num_strips
    strip = Vector{Vertex}[]
    for i in 1:ns
        gvl = unsafe_load(tri_out[].strip, i)
        verts = [unsafe_load(gvl.vertex, k) for k in 1:gvl.num_vertices]
        push!(strip, verts)
    end

    ccall(  (:gpc_free_tristrip, "/usr/local/lib/libgpc.so") , Cvoid,
            (Ptr{gpc_tristrip},), tri_out)

    return [GPCTriStrip(v) for v in strip]

end

"""
`union_strip(p1,p2)`

Computes the union between two [`GPCPolygon`](@ref)s but returns triangle strips instead
of another  [`GPCPolygon`](@ref).
"""
union_strip(p1::GPCPolygon, p2::GPCPolygon) = gpc_tristrip_clip(GPC_UNION, p1, p2)

"""
`intersect_strip(p1,p2)`

Computes the intersection between two [`GPCPolygon`](@ref)s but returns triangle strips
 instead of another  [`GPCPolygon`](@ref).
"""
intersec_strip(p1::GPCPolygon, p2::GPCPolygon) = gpc_tristrip_clip(GPC_INT, p1, p2)


"""
`diff_strip(p1,p2)`

Computes the difference between two [`GPCPolygon`](@ref)s but returns triangle strips 
instead of another  [`GPCPolygon`](@ref).
"""
diff_strip(p1::GPCPolygon, p2::GPCPolygon) = gpc_tristrip_clip(GPC_DIFF, p1, p2)

"""
`xor_strip(p1,p2)`

Computes the exclusive or between two [`GPCPolygon`](@ref)s but returns triangle strips 
instead of another  [`GPCPolygon`](@ref).
"""
xor_strip(p1::GPCPolygon, p2::GPCPolygon) = gpc_tristrip_clip(GPC_XOR, p1, p2)

"""
`pintri(p, p1, p2, p3)`

Find whether point `p` is inside the triangle formed by
vertices `p1`, `p2` and `p3`.
"""
function pintri(p, p1, p2, p3)
    x = p.x-p1.x
    y = p.y-p1.y

    x1 = p2.x-p1.x
    y1 = p2.y-p1.y

    x2 = p3.x-p1.x
    y2 = p3.y-p1.y

    D = x1*y2 - x2*y1
    a =  y2*x - x2*y
    b = -y1*x + x1*y
    if D >= 0
        return a ≥ 0 && b ≥ 0 && (a+b) ≤ D
    else
        return a ≤ 0 && b ≤ 0 && (a+b) ≥ D
    end
    
end


import Base.(∈)
"""
`p ∈ poly`

`p in poly`

`in(p, poly)`

Determine whether point `p` is inside polygon `poly`.

This function will convert the polygon to triangle strips
and check each triangle to see if the point is inside it
using function [`pintri`](@ref).
"""
function (∈)(p::Vertex, poly::GPCPolygon)
    # Convert the polygon to GPCTriStrip
    strip = gpc2tristrip(poly)

    isin = false

    for s in strip
        ntri = length(s)
        for i in 1:ntri
            p1, p2, p3 = s[i]
            if pintri(p, p1, p2, p3)
                return true
            end
        end
    end
    
    return false
end

end
