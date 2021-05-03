#include "Hittable.h"

using namespace std;

// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    glm::vec3 oc = ray.o - o_;
    float a = glm::dot(ray.d, ray.d);
    float b = 2.0 * (glm::dot(oc, ray.d));
    float c = glm::dot(oc, oc) - (r_ * r_);

    float discriminant = (b * b) - (4 * a * c);
    Point intersection(0, 0, 0);
    float t = 0.0;

    if (discriminant > 0) {
        float numerator = -b - glm::sqrt(discriminant);
        if (numerator < 0.0) {
            numerator = -b + glm::sqrt(discriminant);
        }
        t = numerator / (2.0 * a);
        intersection[0] = ray.o[0] + (ray.d[0] * t);
        intersection[1] = ray.o[1] + (ray.d[1] * t);
        intersection[2] = ray.o[2] + (ray.d[2] * t);
        ret = true;
    }

    if (ret) {
        float magnitude = glm::sqrt((intersection[0] * intersection[0]) +
                (intersection[1] * intersection[1]) +
                (intersection[2] * intersection[2]));
        hit_record->normal[0] = intersection[0] / magnitude;
        hit_record->normal[1] = intersection[1] / magnitude;
        hit_record->normal[2] = intersection[2] / magnitude;
        hit_record->in_direction = ray.d;
        // Lecture note 8, P.11
        hit_record->reflection = ((2 * glm::dot(ray.o - intersection, hit_record->normal)) * hit_record->normal) - (ray.o - intersection);
        hit_record->position = intersection;
        hit_record->distance = t;
        hit_record->material = material_;
    }

    return ret;
}

// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    glm::vec4 O = glm::vec4( ray.o, 1.0 );
    glm::vec4 D = glm::vec4( ray.d, 0.0 );
    float a = glm::dot(D, A_ * D);
    float b = 2 * (glm::dot(O, A_ * D));
    float c = glm::dot(O, A_ * O);

    float discriminant = (b * b) - (4 * a * c);
    Point intersection(0, 0, 0);
    float t = 0.0;

    if (discriminant >= 0) {

        if (discriminant == 0) {
            t = -b / 2 * a;
        }
        else {
            float numerator_1 = -b - glm::sqrt(discriminant);
            float numerator_2 = -b + glm::sqrt(discriminant);
            float t_1 = numerator_1 / (2.0 * a);
            float t_2 = numerator_2 / (2.0 * a);
            if (t_1 <= t_2) {
                t = t_1;
            }
            else {
                t = t_2;
            }
        }
        intersection[0] = ray.o[0] + (ray.d[0] * t);
        intersection[1] = ray.o[1] + (ray.d[1] * t);
        intersection[2] = ray.o[2] + (ray.d[2] * t);
        ret = true;
    }

    if (ret) {
        float magnitude = glm::sqrt((intersection[0] * intersection[0]) +
                                    (intersection[1] * intersection[1]) +
                                    (intersection[2] * intersection[2]));
        hit_record->normal[0] = intersection[0] / magnitude;
        hit_record->normal[1] = intersection[1] / magnitude;
        hit_record->normal[2] = intersection[2] / magnitude;
        hit_record->in_direction = ray.d;
        // Lecture note 8, P.11
        hit_record->reflection = ((2 * glm::dot(ray.o - intersection, hit_record->normal)) * hit_record->normal) - (ray.o - intersection);
        hit_record->position = intersection;
        hit_record->distance = t;
        hit_record->material = material_;
    }

    return ret;
}

// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    // TODO: Add your code here.
    bool ret = false;

    glm::vec3 ab = b_-a_;
    glm::vec3 ac = c_-a_;
    glm::vec3 N = glm::cross(ab, ac);

    float n_dot_d = glm::dot(N, ray.d);
    float t = 0.0;
    Point intersection(0, 0, 0);

    if (n_dot_d != 0) {
        float d = glm::dot(N, a_);
        t = (d - glm::dot(N, ray.o)) / n_dot_d;

        intersection[0] = ray.o[0] + (ray.d[0] * t);
        intersection[1] = ray.o[1] + (ray.d[1] * t);
        intersection[2] = ray.o[2] + (ray.d[2] * t);

        // Edge ab
        glm::vec3 edge_ab = b_ - a_;
        glm::vec3 vp_ab = intersection - a_;
        glm::vec3 cross_ab = glm::cross(edge_ab, vp_ab);
        float dot_ab = glm::dot(cross_ab, N);

        // Edge ac
        glm::vec3 edge_bc = c_ - b_;
        glm::vec3 vp_bc = intersection - b_;
        glm::vec3 cross_bc = glm::cross(edge_bc, vp_bc);
        float dot_bc = glm::dot(cross_bc, N);

        // Edge bc
        glm::vec3 edge_ca = a_ - c_;
        glm::vec3 vp_ca = intersection - c_;
        glm::vec3 cross_ca = glm::cross(edge_ca, vp_ca);
        float dot_ca = glm::dot(cross_ca, N);

        if (dot_ab >= 0 && dot_bc >= 0 && dot_ca >= 0) {
            ret = true;
        }
    }

    if (ret) {
        hit_record->in_direction = ray.d;
        // Lecture note 8, P.11
        hit_record->reflection = ((2 * glm::dot(ray.o - intersection, hit_record->normal)) * hit_record->normal) - (ray.o - intersection);
        hit_record->position = intersection;
        hit_record->distance = t;
        if (phong_interpolation_) {
            glm::vec3 ab = b_-a_;
            glm::vec3 ac = c_-a_;
            glm::vec3 N = glm::cross(ab, ac);

            float alpha_1 = glm::dot(glm::cross(c_-b_,intersection-b_), N) /
                    glm::dot(glm::cross(b_-a_,c_-a_), N);
            float alpha_2 = glm::dot(glm::cross(a_-c_,intersection-c_), N) /
                    glm::dot(glm::cross(b_-a_,c_-a_), N);
            float alpha_3 = glm::dot(glm::cross(b_-a_,intersection-a_), N) /
                    glm::dot(glm::cross(b_-a_,c_-a_), N);
            glm::vec3 temp_normal = (alpha_1 * n_a_) + (alpha_2 * n_b_) + (alpha_3 * n_c_);
            hit_record->normal = glm::normalize(temp_normal);

        }
        else {
            glm::vec3 edge_ab = b_ - a_;
            glm::vec3 edge_ac = c_ - a_;
            hit_record->normal = glm::cross(edge_ab, edge_ac);

        }
        // no need to set material in this function
    }
    return ret;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord *hit_record) const {
    bool ret = triangle_.Hit(ray, hit_record);
    if (ret) {
        hit_record->material = material_;
    }
    return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
           const Material& material,
           bool phong_interpolation):
           ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
    std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
    vertices_.resize(v_pos.size());

    for (int i = 0; i < vertices_.size(); i++) {
        vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
    }

    f_ind_ = ply_data_.getFaceIndices();

    // Calc face normals
    for (const auto& face : f_ind_) {
        Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
        face_normals_.emplace_back(normal);
    }

    // Calc vertex normals
    vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
    for (int i = 0; i < f_ind_.size(); i++) {
        for (int j = 0; j < 3; j++) {
            vertex_normals_[f_ind_[i][j]] += face_normals_[i];
        }
    }
    for (auto& vertex_normal : vertex_normals_) {
        vertex_normal = glm::normalize(vertex_normal);
    }

    // Construct hittable triangles
    for (const auto& face : f_ind_) {
        triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
                                vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
                                phong_interpolation_);
    }

    // Calc bounding box
    Point bbox_min( 1e5f,  1e5f,  1e5f);
    Point bbox_max(-1e5f, -1e5f, -1e5f);
    for (const auto& vertex : vertices_) {
        bbox_min = glm::min(bbox_min, vertex - 1e-3f);
        bbox_max = glm::max(bbox_max, vertex + 1e-3f);
    }

    // Build Octree
    tree_nodes_.emplace_back(new OctreeNode());
    tree_nodes_.front()->bbox_min = bbox_min;
    tree_nodes_.front()->bbox_max = bbox_max;

    root_ = tree_nodes_.front().get();
    for (int i = 0; i < f_ind_.size(); i++) {
        InsertFace(root_, i);
    }
}

bool Mesh::Hit(const Ray& ray, HitRecord *hit_record) const {
    const bool brute_force = false;
    if (brute_force) {
        // Naive hit algorithm
        float min_dist = 1e5f;
        for (const auto &triangle : triangles_) {
            HitRecord curr_hit_record;
            if (triangle.Hit(ray, &curr_hit_record)) {
                if (curr_hit_record.distance < min_dist) {
                    *hit_record = curr_hit_record;
                    min_dist = curr_hit_record.distance;
                }
            }
        }
        if (min_dist + 1.0 < 1e5f) {
            hit_record->material = material_;
            return true;
        }
        return false;
    } else {
        bool ret = OctreeHit(root_, ray, hit_record);
        if (ret) {
            hit_record->material = material_;
        }
        return ret;
    }
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
    for (size_t idx : face) {
        const auto& pt = vertices_[idx];
        for (int i = 0; i < 3; i++) {
            if (pt[i] < bbox_min[i] + 1e-6f) return false;
            if (pt[i] > bbox_max[i] - 1e-6f) return false;
        }
    }
    return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
    float t_min = -1e5f;
    float t_max =  1e5f;

    for (int i = 0; i < 3; i++) {
        if (glm::abs(ray.d[i]) < 1e-6f) {
            if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
                t_min =  1e5f;
                t_max = -1e5f;
            }
        }
        else {
            if (ray.d[i] > 0.f) {
                t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
                t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
                t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
            }
        }
    }

    return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
    const Point& bbox_min = u->bbox_min;
    const Point& bbox_max = u->bbox_max;

    Vec bias = bbox_max - bbox_min;
    Vec half_bias = bias * 0.5f;

    bool inside_childs = false;

    for (size_t a = 0; a < 2; a++) {
        for (size_t b = 0; b < 2; b++) {
            for (size_t c = 0; c < 2; c++) {
                size_t child_idx = ((a << 2) | (b << 1) | c);
                Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
                Point curr_bbox_max = curr_bbox_min + half_bias;
                if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
                    if (u->childs[child_idx] == nullptr) {
                        tree_nodes_.emplace_back(new OctreeNode());
                        OctreeNode* child = tree_nodes_.back().get();
                        u->childs[child_idx] = tree_nodes_.back().get();
                        child->bbox_min = curr_bbox_min;
                        child->bbox_max = curr_bbox_max;
                    }
                    InsertFace(u->childs[child_idx], face_idx);
                    inside_childs = true;
                }
            }
        }
    }

    if (!inside_childs) {
        u->face_index.push_back(face_idx);
    }
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
    if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
        return false;
    }
    float distance = 1e5f;
    for (const auto& face_idx : u->face_index) {
        HitRecord curr_hit_record;
        if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }

    for (const auto& child : u->childs) {
        if (child == nullptr) {
            continue;
        }
        HitRecord curr_hit_record;
        if (OctreeHit(child, ray, &curr_hit_record)) {
            if (curr_hit_record.distance < distance) {
                distance = curr_hit_record.distance;
                *hit_record = curr_hit_record;
            }
        }
    }
    return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
    hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord *hit_record) const {
    float min_dist = 1e5f;
    for (const auto &hittable : hittable_list_) {
        HitRecord curr_hit_record;
        if (hittable->Hit(ray, &curr_hit_record)) {
            if (curr_hit_record.distance < min_dist) {
                *hit_record = curr_hit_record;
                min_dist = curr_hit_record.distance;
            }
        }
    }
    return min_dist + 1.0 < 1e4f;
}