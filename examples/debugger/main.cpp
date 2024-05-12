#include "Delaunay_psm.h"
#include "mantis.h"
#include "util.h"

#include <random>

namespace ps = polyscope;

ps::SurfaceMesh *tool = nullptr;
ps::SurfaceMesh *mesh = nullptr;

// quad face in xz plane
std::vector<glm::vec3> points = {{-10, 0, 10},
                                 {10,  0, 10},
                                 {10,  0, -10},
                                 {-10, 0, -10}};
std::vector<glm::vec3> points_og = {{-10, 0, 10}, {10, 0, 10}, {10, 0, -10}, {-10, 0, -10}};
std::vector<std::array<int, 4>> faces = {{0, 1, 2, 3}};

class DisjointSets
{
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    explicit DisjointSets(int size)
    {
        parent.resize(size);
        rank.resize(size, 0);

        for (int i = 0; i < size; ++i)
        {
            parent[i] = i;
        }
    }

    int find(int x)
    {
        if (parent[x] != x)
        {
            parent[x] = find(parent[x]); // Path compression
        }
        return parent[x];
    }

    void unionSets(int x, int y)
    {
        int xRoot = find(x);
        int yRoot = find(y);

        if (xRoot == yRoot)
        {
            return;
        }

        // Union by rank
        if (rank[xRoot] < rank[yRoot])
        {
            parent[xRoot] = yRoot;
        }
        else if (rank[yRoot] < rank[xRoot])
        {
            parent[yRoot] = xRoot;
        }
        else
        {
            parent[yRoot] = xRoot;
            rank[xRoot] = rank[xRoot] + 1;
        }
    }
};

double eval_plane(GEO::vec4 plane, GEO::vec3 p) {
    return plane.x * p.x + plane.y * p.y + plane.z * p.z + plane.w;
}

struct ConvexCellV1 {

    ConvexCellV1() = default;

    ConvexCellV1(const std::vector<glm::vec3>& pts, const std::vector<std::array<int, 4>>& quads) {
        points.resize(pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            points[i] = {pts[i].x, pts[i].y, pts[i].z};
        }

        std::map<std::pair<index_t, index_t>, std::pair<index_t, index_t>> edge_map;

        for (size_t f = 0; f < quads.size(); ++f) {
            auto q = quads[f];
            for (size_t i = 0; i < 4; ++i) {
                auto key = std::minmax(q[i], q[(i + 1) % 4]);
                auto& p = edge_map[key];
                if(key.first == q[i]) {
                    p.first = f;
                } else {
                    p.second = f;
                }
            }
        }

        for(auto e: edge_map) {
            edges.push_back(ConvexCellV1::Edge{e.first.first, e.first.second, e.second.first, e.second.second});
        }

        original_edges = edges;
        num_faces = quads.size();
        original_num_faces = num_faces;
        original_num_vertices = points.size();
    }

    struct Edge {
        index_t from = -1, to = -1;
        index_t left = -1, right = -1;
    };

    std::vector<GEO::vec3> points;

    std::vector<Edge> edges;
    std::vector<Edge> original_edges;

    // temporary buffers
    mutable std::vector<Edge> edge_buffer;
    mutable std::vector<int> vertices_buffer;
    mutable std::vector<double> distance;
    mutable std::vector<int> label;

    index_t num_faces = 0;
    index_t original_num_faces = 0;
    index_t original_num_vertices = 0;

    void clip(GEO::vec4 plane) {
        edge_buffer.clear();
        vertices_buffer.clear();

        distance.resize(points.size());
        label.resize(points.size());

        bool has_positive = false;
        bool has_negative = false;

        //vertices.resize(points.size());
        //std::iota(vertices.begin(), vertices.end(), 0);

        for (int v = 0; v < points.size(); ++v) {
            double d = eval_plane(plane, points[v]);
            has_positive |= d > 0;
            has_negative |= d <= 0;
            distance[v] = d;
            label[v] = d > 0 ? 1 : -1;
        }

        if (!has_negative) {
            //printf("no negative\n");
            return;
        }

        if (!has_positive) {
            // result is empty
            edges.clear();
            return;
        }

        //make_consistent();

        for(int v = 0; v < points.size(); ++v) {
            if(label[v] > 0) {
                vertices_buffer.push_back(v);
            }
        }

        std::vector<index_t> x_edge(num_faces, -1);
        index_t new_face_idx = num_faces++;

        auto try_make_edge = [&](index_t f, index_t v, bool is_from) {
            if (x_edge[f] == -1) {
                x_edge[f] = v;
            } else {
                Edge e{.left = new_face_idx, .right = f, .from = is_from ? x_edge[f] : v, .to = is_from ? v
                                                                                                        : x_edge[f]};
                edge_buffer.push_back(e);
            }
        };

        for (Edge& edge : edges) {
            index_t from = edge.from;
            index_t to = edge.to;

            double from_dist = distance[from];
            double to_dist = distance[to];

            int from_label = label[from];
            int to_label = label[to];

            if (from_label * to_label < 0) {
                double t = from_dist / (from_dist - to_dist);
                GEO::vec3 pt = points[from] + t * (points[to] - points[from]);
                auto v = points.size();
                points.push_back(pt);
                bool is_from = from_label > 0;
                if (from_label > 0) {
                    edge.to = v;
                } else {
                    edge.from = v;
                }
                edge_buffer.push_back(edge);
                try_make_edge(edge.left, v, is_from);
                try_make_edge(edge.right, v, !is_from);
            }
            else if(from_label > 0 && to_label > 0) {
                edge_buffer.push_back(edge);
            }
        }

        std::swap(edges, edge_buffer);
    }

    void make_consistent() {
        int num_vertices = points.size();
        DisjointSets ds(num_vertices);
        for (const Edge& edge : edges) {
            if(label[edge.from] == label[edge.to]) {
                ds.unionSets(edge.from, edge.to);
            }
        }

        std::vector<int> component_size(num_vertices, 0);
        int num_components = 0;
        for(int i = 0; i < num_vertices; ++i) {
            int root = ds.find(i);
            if(root == i) {
                num_components++;
            }
            component_size[root]++;
        }

        if(num_components == 2) {
            return;
        }

        int largest_component_a = -1, largest_component_b = -1;
        for(int i = 0; i < num_vertices; ++i) {
            if(component_size[i] > 0) {
                if(label[i] > 0 && (largest_component_a == -1 || component_size[i] > component_size[largest_component_a])) {
                    largest_component_a = i;
                } else if(label[i] <= 0 && (largest_component_b == -1 || component_size[i] > component_size[largest_component_b])) {
                    largest_component_b = i;
                }
            }
        }

        if(largest_component_a == -1 || largest_component_b == -1) return; // No need to proceed if we can't find two components.

        std::vector<int> starts(num_vertices + 1, 0);
        for(const Edge& edge : edges) {
            starts[edge.from + 1]++;
            starts[edge.to + 1]++;
        }
        std::partial_sum(starts.begin(), starts.end(), starts.begin());
        std::vector<int> neighbors(starts.back()), position = starts;

        for(const Edge& e : edges) {
            neighbors[position[e.from]++] = e.to;
            neighbors[position[e.to]++] = e.from;
        }

        auto flood_fill = [&](int start_vertex, int avoid_component) {
            std::vector<bool> visited(num_vertices, false);
            std::vector<int> queue = {start_vertex};
            visited[start_vertex] = true;

            while (!queue.empty()) {
                int current = queue.back();
                queue.pop_back();
                for (int i = starts[current]; i < starts[current + 1]; ++i) {
                    int neighbor = neighbors[i];
                    if (!visited[neighbor] && ds.find(neighbor) != avoid_component) {
                        visited[neighbor] = true;
                        ds.unionSets(current, neighbor);
                        queue.push_back(neighbor);
                    }
                }
            }
        };

        flood_fill(largest_component_a, largest_component_b);
        flood_fill(largest_component_b, largest_component_a);
    }

    void reset() {
        edges = original_edges;
        num_faces = original_num_faces;
        points.resize(original_num_vertices);
    }
};

struct ConvexCell {
    struct Halfedge {
        int next;
        int face;
        int vertex; // tip vertex
    };

    std::vector<Halfedge> halfedges;
    std::vector<int> vertex_halfedge; // incoming he
    std::vector<int> face_halfedge;
    std::vector<GEO::vec3> points;

    void clip(GEO::vec4 plane) {
        size_t nv = (int)vertex_halfedge.size();

        //std::vector<int> component(nv, -1);

        // find connected components by traversing the mesh using the halfedge connectivity
        //std::vector<int> q;
        //size_t num_vertices_visited = 0;
        //int num_components = 0;
        //auto traverse_from_vertex = [&](int start_vertex) {
        //    q.clear();
        //    q.push_back(int(start_vertex));
        //    int c = num_components++;
        //    while(!q.empty()) {
        //        int v = q.back();
        //        q.pop_back();
        //        if(component[v] != -1) continue;

        //        // assign component id
        //        ++num_vertices_visited;
        //        component[v] = c;

        //        // iterate over outgoing halfedges
        //        int he = vertex_halfedge[v];
        //        int start = he;
        //        do {
        //            int next = halfedges[he].next;
        //            int to = halfedges[he].vertex;
        //            if(label[to] == label[v]) {
        //                q.push_back(to);
        //            }
        //            else {
        //                mixed_he = he;
        //                //mixed_edges.push_back(he);
        //            }
        //            he = next ^ 1;
        //        } while(he != start);
        //    }
        //};

        //int v = traverse_from_vertex(0);
        //traverse_from_vertex(v);

        //if(num_vertices_visited == nv) {
        //    // do the cut

        //}

        //if(num_components > 2) {
        //    ensure_two_connected_components(component, major_pos, major_neg);
        //}

        auto vertex_dist = [&](int v) {
            return eval_plane(plane, points[v]);
        };

        auto vertex_component = [&](int v) {
            return vertex_dist(v) > 0 ? 1 : -1;
        };

        // find a mixed edge
        int mixed_he = -1;
        for(int e = 0; e < halfedges.size(); ++e) {
            int from = halfedges[e].vertex;
            int to = halfedges[e ^ 1].vertex;
            int from_c = vertex_component(from);
            int to_c = vertex_component(to);
            if(from_c * to_c < 0) {
                mixed_he = e;
                break;
            }
        }

        int nf = (int)face_halfedge.size();

        std::vector<index_t> x_edge(nf, -1);
        int new_face_idx = nf;

        //int keep = component[major_pos];

        // orient mixed halfedge outwards, so the halfedge vertex is negative
        if(vertex_dist(halfedges[mixed_he].vertex) < 0 ){
            mixed_he = mixed_he ^ 1;
        }
        int he = mixed_he;

        auto find_next_mixed_edge = [&](int he){
            while(true){
                he = halfedges[he].next;
                //int c_from = component[halfedges[he].vertex];
                int c_from = vertex_component(halfedges[he].vertex);
                //int c_to = component[halfedges[he ^ 1].vertex];
                int c_to = vertex_component(halfedges[he ^ 1].vertex);
                if(c_from != c_to) {
                    break;
                }
            }
            return he ^ 1;
        };

        do {
            int from = halfedges[he ^ 1].vertex;
            int to = halfedges[he].vertex;

            double from_dist = vertex_dist(from);
            double to_dist = vertex_dist(to);

            double t = from_dist / (from_dist - to_dist);
            GEO::vec3 pt = points[from] + t * (points[to] - points[from]);

            int new_vertex = (int)points.size();
            points.push_back(pt);
            int next_mixed = find_next_mixed_edge(he);

            int f = halfedges[he].face;
            halfedges[he].vertex = new_vertex;
            face_halfedge[f] = he;
            vertex_halfedge[new_vertex] = he;
            halfedges.push_back(Halfedge{.next = next_mixed ^ 1, .face = f, .vertex = new_vertex});
            int last_outer_he = halfedges.size() - 2;
            halfedges.push_back(Halfedge{.next = last_outer_he, .face = new_face_idx, .vertex = to});

            he = next_mixed;
        } while(he != mixed_he);
    }

    void ensure_two_connected_components(std::vector<int>& component, int major_a, int major_b) {
        size_t nv = points.size();
        auto flood_fill = [&](int start_vertex, int avoid_component) {
            std::vector<bool> visited(nv, false);
            std::vector<int> q = {start_vertex};
            visited[start_vertex] = true;

            while (!q.empty()) {
                int v = q.back();
                q.pop_back();

                // iterate over outgoing halfedges
                int he = vertex_halfedge[v];
                int start = he;
                do {
                    int next = halfedges[he].next;
                    int to = halfedges[he].vertex;
                    if(component[to] != avoid_component) {
                        q.push_back(to);
                    }
                    he = next ^ 1;
                } while(he != start);
            }
        };

        flood_fill(major_a, component[major_b]);
        flood_fill(major_b, component[major_a]);
    }
};

ConvexCellV1 cell;

void draw_cellv1(const ConvexCellV1& cell) {
    std::vector<glm::vec3> pts;

    pts.reserve(cell.points.size());
    for(auto p: cell.points) {
        pts.emplace_back(p.x, p.y, p.z);
    }

    //compute face_centroids
    std::map <index_t, std::pair<glm::vec3, int>> face_centroid;

    //int count = 0;
    // for each face find some vertex
    for (auto e: cell.edges) {
        auto v0 = e.from;
        auto v1 = e.to;
        //draw_edge("edge" + std::to_string(count++), cell.points[v0], cell.points[v1]);
        {
            auto [c, n] = face_centroid[e.left];
            face_centroid[e.left] = {c + pts[v0] + pts[v1], n + 2};
        }

        {
            auto [c, n] = face_centroid[e.right];
            face_centroid[e.right] = {c + pts[v0] + pts[v1], n + 2};
        }
    }

    std::map<index_t, index_t> face_vertex;
    // use face centroid as face vertex
    for(auto [f, p]: face_centroid) {
        auto [c, n] = p;
        c /= n;
        pts.push_back(c);
        face_vertex[f] = pts.size() - 1;
    }

    std::vector<std::array<index_t, 3>> tris;

    // generate fan triangulation for each face
    for (auto e: cell.edges) {
        auto v0 = e.from;
        auto v1 = e.to;
        {
            assert(face_vertex.count(e.left) > 0);
            auto v2 = face_vertex[e.left];
            if(v2 != v0 && v2 != v1) {
                tris.push_back({v0, v1, v2});
            }
        }

        {
            assert(face_vertex.count(e.right) > 0);
            auto v2 = face_vertex[e.right];
            if(v2 != v0 && v2 != v1) {
                tris.push_back({v1, v0, v2});
            }
        }
    }

    if(tris.empty()) {
        pts.clear();
    }

    ps::registerSurfaceMesh("cell", pts, tris);
}

// for a vertex-element combination render the interception region
void callback() {
    glm::mat4 tf = tool->getTransform();

    if (ImGui::Button("Do Cut")) {
        // xz plane
        glm::vec4 plane(0, 1, 0, 0);
        // transform with tf
        plane = glm::transpose(glm::inverse(tf)) * plane;
        cell.clip({plane.x, plane.y, plane.z, plane.w});
        draw_cellv1(cell);
    }

    if(ImGui::Button("Reset")) {
        cell.reset();
        draw_cellv1(cell);
    }
}

void load_box(std::vector<glm::vec3> &pts, std::vector<std::array<int, 4>> &quads) {

    glm::vec3 lower = {-1, -1, -1};
    glm::vec3 upper = {1, 1, 1};

    pts
            = { /* lower face */ lower,        { upper.x, lower.y, lower.z },
                                 { upper.x, upper.y, lower.z }, { lower.x, upper.y, lower.z },
                    /* upper face */ upper,        { lower.x, upper.y, upper.z },
                                 { lower.x, lower.y, upper.z }, { upper.x, lower.y, upper.z } };

    quads.push_back({ 0, 3, 2, 1 });
    quads.push_back({ 4, 5, 6, 7 });
    quads.push_back({ 0, 6, 5, 3 });
    quads.push_back({ 1, 2, 4, 7 });
    quads.push_back({ 1, 7, 6, 0 });
    quads.push_back({ 2, 3, 5, 4 });
}

int main(int, char **) {
    ps::init();
    ps::options::giveFocusOnShow = true;

    std::vector<glm::vec3> pts;
    std::vector<std::array<int, 4>> quads;
    load_box(pts, quads);

    //glm::mat3 random_rot = glm::mat3_cast(glm::angleAxis(0.5f, glm::vec3(0.2, 0.3, 0.7)));
    //for(auto &p: pts) {
    //    p = random_rot * p;
    //}

    // transform plane as well
    //for(auto &p: points) {
    //    p = random_rot * p;
    //}


    cell = ConvexCellV1(pts, quads);

    draw_cellv1(cell);

    tool = ps::registerSurfaceMesh("tool", points, faces);
    tool->setBackFaceColor({1, 0, 1});
    tool->setBackFacePolicy(ps::BackFacePolicy::Custom);

    ps::state::userCallback = callback;
    ps::show();
}
