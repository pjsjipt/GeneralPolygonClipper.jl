module GeneralPolygonClipper
import GeometryBasics: Point, Point2, Triangle, Polygon, AbstractPolygon


@enum GPCOperation GPC_DIFF=0 GPC_INT GPC_XOR GPC_UNION


export Point2, Polygon, Triangle, area
export GPCPolygon, gpc_polygon_clip, gpc_tristrip_clip
export union_strip, intersect_strip, diff_strip, xor_strip
export gpcpoly2tristrip, tristrip2triangles, tristrip2polygons
export gpcpoly2triangles, gpcpoly2polygons
export gpc2poly
export GPCOperation, GPC_DIFF, GPC_INT, GPC_XOR, GPC_UNION

const Vertex = Point2{Float64}

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
struct GPCPolygon <: AbstractPolygon{2,Float64}
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


import GeometryBasics.area

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

`gpcpoly2tristrip(p::GPCPolygon)`

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

    return strip
    
end


"""
`gpc_polygon_to_tristrip(p::GPCPolygon)`

`gpcpoly2tristrip(p::GPCPolygon)`

Converts a [`GPCPolygon`](@ref) to a triangle strip.

In GPC, a triangle strip is a vector of vectors of vertices that specify a triangle strip. Actually, each vector of vertices correspond to a contiguous triangle strip and a polygon can be decomposed into several triangle strip. 

Thus this function will return a vector of vector of vertices.

In each vector of vertices, the 3 first vertices make up a triangle. Each additional
vertex adds a new triangle to the strip. See function [`trianglevertices`](@ref) to see how each triangle is made up

"""
gpcpoly2tristrip(p::GPCPolygon) = gpc_polygon_to_tristrip(p)

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

"""
`tristrip2triangles(strip)`

Convert each individual triangle strip to a vector of [`GeometryBasics::Triangle`](@ref) objects in [`GeometryBasics`](@ref) package.
"""
function tristrip2triangles(strip::Vector{P}) where {P <: Point}

    nv = length(strip)
    @assert nv >= 3
    fun(idx) = Triangle(strip[idx[1]], strip[idx[2]], strip[idx[3]])
    return [fun(trianglevertices(ivert)) for ivert in 1:nv-2]
end

"""
`tristrip2polygons(strip)`

Convert each individual triangle strip to a vector of [`GeometryBasics::Polygon`](@ref) objects in [`GeometryBasics`](@ref) package.
"""
function tristrip2polygons(strip::Vector{P}) where {P <: Point}
    nv = length(strip)
    @assert nv >= 3
    fun(idx) = Polygon([strip[idx[1]], strip[idx[2]], strip[idx[3]]])
    return [fun(trianglevertices(ivert)) for ivert in 1:nv-2]
end

    
"""
`gpcpoly2triangles(p)`

Decompose a 
"""
function gpcpoly2triangles(p::GPCPolygon)
    strips = gpcpoly2tristrip(p)

    nstrips = length(strips)

    if nstrips == 0
        return Triangle{2,Float64}[]
    elseif nstrips == 1
        return tristrip2triangles(strips[1])
    else
        return reduce(append!, [tristrip2triangles(s) for s in strips])
    end
end

function gpcpoly2polygons(p::GPCPolygon)
    strips = gpcpoly2tristrip(p)

    nstrips = length(strips)

    if nstrips == 0
        return typeof(Polygon([Vertex(0.0,0.0)]))
    elseif nstrips == 1
        return tristrip2polygons(strips[1])
    else
        return reduce(append!, [tristrip2triangles(s) for s in strips])
    end
end


"""
`gpc2poly(p::GPCPolygon)` 

Converts a [`GPCPolygon`](@ref) object to a 
[`GeometryBasics::Polygon`](@ref).

The input polygon should have one and only one external contour!
"""
function gpc2poly(p::GPCPolygon)
    # This will only work *if* we have a single external contour
    if sum(!, p.holes) != 1
        error("GPCPolygon should have only one external contour!")
    end
    
    nc = length(p.contours)
    cext = [p.contours[i] for i in 1:nc if !p.holes[i]]
    cint = [p.contours[i] for i in 1:nc if p.holes[i]]

    if length(cint) == 0
        return Polygon(cext[1])
    else
        return Polygon(cext[1], cint)
    end
    
end



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

    return strip

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


end
