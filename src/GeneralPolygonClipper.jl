module GeneralPolygonClipper
import GeometryBasics: Point2

const Vertex = Point2{Float64}


struct gpc_vertex_list
    num_vertices::Cint
    vertex::Ptr{Vertex}
end



struct gpc_polygon
    num_contours::Cint
    hole::Ptr{Cint}
    contour::Ptr{gpc_vertex_list}
end

struct GPCPolygon
    holes::Vector{Bool}
    contours::Vector{Vector{Vertex}}
end

struct gpc_tristrip
    num_strips::Cint
    strip::Ptr{gpc_vertex_list}
end

function GPCPolygon(x::AbstractVector, y::AbstractVector)
    nv = length(x)
    @assert nv == length(y)

    v = [Vertex(x[i], y[i]) for i in 1:nv]

    return GPCPolygon([false], [v])
end


#const plibgpc = Libc.Libdl.dlopen("/usr/local/lib/libgpc.so")
#const pclip = Libc.Libdl.dlsym(plibgpc, "gpc_polygon_clip")
#const pfree = Libc.Libdl.dlsym(plibgpc, "gpc_free_polygon")
#const pteste = Libc.Libdl.dlsym(plibgpc, "teste")
#const plixo = Libc.Libdl.dlsym(plibgpc, "caralho")
teste(x) = ccall( (:teste, "/usr/local/lib/libgpc.so") , Cint, (Cint,), x)


function gpc_polygon_clip(op::Integer, p1::GPCPolygon, p2::GPCPolygon)

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


function gpc_polygon_to_tristrip(p::GPCPolygon)
    nc1 = Cint(length(p.holes))
    h1 = [Cint(h) for h in p.holes]
    c1 = [gpc_vertex_list(length(vlst), pointer(vlst)) for vlst in p.contours]
    gpc1 = gpc_polygon(nc1, pointer(h1), pointer(c1))

    tri_out = Ref(gpc_tristrip(0, 0))

    ccall( (:gpc_polygon_to_tristrip, "/usr/local/lib/libgpc.so"), Cvoid,
           (Ptr{gpc_polygon}, Ptr{gpc_tristrip}), Ref(gpc1), tri_out)

    ns = tri_out[].num_strips
    println(ns)
    strip = Vector{Vertex}[]
    for i in 1:ns
        gvl = unsafe_load(tri_out[].strip, i)
        println(gvl.num_vertices)
        verts = [unsafe_load(gvl.vertex, k) for k in 1:gvl.num_vertices]
        push!(strip, verts)
    end

    ccall(  (:gpc_free_tristrip, "/usr/local/lib/libgpc.so") , Cvoid,
            (Ptr{gpc_tristrip},), tri_out)

    return strip
    
end

end
