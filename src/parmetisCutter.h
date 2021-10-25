#ifndef PARMETISDECOMPOSER_H
#define PARMETISDECOMPOSER_H

#include "meshCutter.h"
#include "CGFile.h"

#include <cgnslib.h>
#include <parmetis.h>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <memory>

namespace MeshCut
{
    class ParMetisMeshCutter : public MeshCutter
    {
    public:
        ParMetisMeshCutter(int argc, char** argv);
        virtual ~ParMetisMeshCutter();
    public:
        void cut(const std::string &mesh, const int npart) override;
    
    private:
        using cgt = cgsize_t;
        using vcgt = std::vector<cgt>;
        using vvcgt = std::vector<vcgt>;
        using FileName = std::string;

    public:
        // 
        struct Face_
        {
            // face type
            // cur version only supports below 3 type:
            // TRI_3
            // QUAD_4
            // MIX
            // ElementType_t type = ElementType_t::ElementTypeNull;
            // cgsize_t id = 0;
            // cgsize_t nbrPart = 0;
            std::set<cgsize_t> nodeList;

            bool operator==(const Face_ &other) const
            {
                if (this->nodeList.empty() || other.nodeList.empty())
                    return false;

                auto len = std::min(nodeList.size(), other.nodeList.size());
                auto it = nodeList.begin();
                auto jt = other.nodeList.begin();
                for (auto i = 0; i < len; ++i, ++it, ++jt)
                {
                    if (*it != *jt)
                        return false;
                }
                return true;
            }
        };

        struct FaceInfo
        {
            // face type
            // cur version only supports below 3 type:
            // TRI_3
            // QUAD_4
            // MIX
            ElementType_t type = ElementType_t::ElementTypeNull;

            // face size of each part
            cgsize_t size_all_rank = 0;
            cgsize_t size_this_rank= 0;

            // face id of each part
            cgsize_t id_start_this_rank = 0;
            cgsize_t id_end_this_rank = 0;

            // node count include repeated 
            cgsize_t node_size_this_rank = 0;
            std::vector<cgsize_t> node_id_offset;
        };

    private:
        // 
        struct SectionPartition;

        struct NodeInfo
        {
            // node size of each part
            cgsize_t size_all_rank = 0;
            cgsize_t size_this_rank= 0;

            // node id of each part
            cgsize_t id_start_this_rank = 0;
            cgsize_t id_end_this_rank = 0;

            // original mesh node id range
            cgsize_t id_start_global = 0;
            cgsize_t id_end_global = 0;
        };

        struct CellInfo
        {
            // cell size of each part
            cgsize_t size_all_rank = 0;
            cgsize_t size_this_rank= 0;

            // cell id of each part
            cgsize_t id_start_this_rank = 0;
            cgsize_t id_end_this_rank = 0;
        };

    private:
        // parallel rank/size 
        int size_=1, rank_=0;
        bool isMaster_ = false;

        // decompose number
        idx_t nparts_=1;
        std::shared_ptr<CGIO::CGFile> bigFile_;
        std::vector<std::shared_ptr<CGIO::CGFile>> smallFile_;

        // section 
        std::vector<cgt> bodySections_;
        std::vector<cgt> boundarySections_;
        std::vector<SectionPartition> sectionPartition;
        
        // cell/node info 
        cgt totalCell_=0, totalNode_=0;
        std::vector<std::vector<NodeInfo>> nodeInfo_;
        std::vector<FaceInfo> faceInfo_;
        std::vector<std::vector<CellInfo>> cellInfo_;
        std::vector<std::unordered_set<Face_>> outsideFaceCell_;
        std::vector<std::unordered_map<cgsize_t, cgsize_t>> outsideNode_;
        std::vector<std::vector<std::vector<std::vector<Face_>>>> interface_;
        std::vector<std::vector<std::unordered_map<cgsize_t, cgsize_t>>> cellPartIdMapToGlobalId_;
        std::vector<std::vector<std::map<cgsize_t, cgsize_t>>> nodeGlobalIdMapToPartId_, nodePartIdMapToGlobalId_;

    private:
        void decompose_body_section();
        void rebuild_cgnsFile(const std::string meshName);
        inline void get_zone_cell_node(const std::vector<cgsize_t>& sectionIDs);
        inline void get_boundary_cell(std::vector<cgsize_t>&& sectionIDs, const cgsize_t ithPart);
        inline std::set<std::set<cgsize_t>> get_outside_faceCell(const cgsize_t ithSection, const cgsize_t ithPart);

        inline void load_coordinate(const FileName& meshName, std::vector<std::vector<std::vector<double>>>& localCoords);
        inline void find_interface_perThread();
        inline void map_interface_allThread();
        inline void map_node_id();
        void charArrayBitOR(char *in, char *inout, int *len, MPI_Datatype *dptr);
        bool checkCGNSFile(const std::shared_ptr<CGIO::CGFile> fp);
    };
    
} // namespace MeshCut

template <>
struct std::hash<MeshCut::ParMetisMeshCutter::Face_>
{
    size_t operator()(const MeshCut::ParMetisMeshCutter::Face_ &k) const
    {
        int h = 1;
        for (auto x : k.nodeList)
        {
            h ^= x;
        }
        return h;
    }
};

#endif
