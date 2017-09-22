#ifndef Triangle_H
#define Triangle_H

#include "ia/edge.h"

class Triangle
{
public:
    Triangle()
    {
        this->vertices = { -1, -1, -1 };
        this->adj = { -1, -1, -1 };
    }

    Triangle(itype v1, itype v2, itype v3)
    {
        this->vertices = { v1, v2, v3 };
        this->adj = { -1, -1, -1 };
    }

    inline itype TV(int pos) { return this->vertices[pos]; }

    inline Edge TE(int pos) { return Edge(vertices[(pos+1)%3],vertices[(pos+2)%3]); }

    inline itype TT(int pos) { return this->adj[pos]; }

    inline void set_TT(int pos, itype adjId) { this->adj[pos]=adjId; }

    inline int vertices_num() { return 3; }

    inline bool operator== (const Triangle &p)
    {
        bool b[3];
        b[0] = false; b[1] = false; b[2] = false;
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(!b[j] && p.vertices[i]==this->vertices[j])
                {
                    b[j] = true;
                    break;
                }
            }
        }

        return b[0] && b[1] && b[2];
    }

    inline bool operator!= (const Triangle p)
    {
       return !((*this)==(p));
    }

    //return the index of vertex v in triangle t
    inline itype vertex_index(itype v)
    {
        for(int i=0; i<3; i++)
            if (TV(i)==v)
                return i;
        return -1;
    }

    inline bool has_vertex(itype v) { return (vertex_index(v)!=-1); }

    //return the next triangle incident in vertex v, starting from triangle 'current'
    inline void next_triangle_around_v(itype center, itype &current, itype &pred)
    {
        int c_pos = this->vertex_index(center);
        if(this->TT((c_pos+1)%3) == pred)
        {
            pred = current;
            current = this->TT((c_pos+2)%3);
        }
        else if(this->TT((c_pos+2)%3) == pred)
        {
            pred = current;
            current = this->TT((c_pos+1)%3);
        }
    }

    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Triangle& p)
    {
        out <<"TV[" << p.vertices[0] << " " << p.vertices[1] << " " << p.vertices[2] << "] ";
        out <<"TT[" << p.adj[0] << " " << p.adj[1] << " " << p.adj[2] << "] ";
        return out;
    }

private:
    ivect vertices;
    ///
    ivect adj;
};



#endif // Triangle_H
