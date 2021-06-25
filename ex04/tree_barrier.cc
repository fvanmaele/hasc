#include "ctb.hh"
#include <cassert>

class CTB
{
private:
    std::condition_variable cv_active;
    //std::condition_variable cv_passive;
    std::mutex mx;
    std::vector<int> flag_grp; // labels (synced groups)
    int P;
    int kn;
    bool tree_ready = false;
    std::vector<int> queue;

public:
    CTB(int P_) : 
        P(P_) 
    {
        kn = P / 2; // number of thread groups
        assert(P == kn*2); // fixme: uneven number of threads
        flag_grp = std::vector<int>(kn, -1); // leaf nodes
    }
    void wait(int rank)
    {
        // xxx: allocated for each thread (stack seperate)
        // -> needs one node per neighboring ranks
        // use class array of cv/mutex instead? (exchange data between threads)
        // but: is this "local" to (groups of) threads? (NUMA-aware)
        auto node = CTBNode(this->tree_ready);
        int i, j;

        if (rank % 2 == 0) {
            i = rank;
            j = rank+1;
        } else {
            i = rank+1;
            j = rank;
        }
        if (j >= P)
            // fixme: case with uneven threads
            std::terminate();

        // check which group we are in
        int k_cur = (rank+1) / kn;

        // returns when i, j have synchronized
        // up: arrived last of (i, j) -> move up the tree
        int up = node.wait(rank, i, j);
        std::unique_lock<std::mutex> ul{ mx }; // mutex on flag_grp
        flag_grp[k_cur] = up;

        int k_nbr;
        if (k_cur % 2 == 0)
            k_nbr = k_cur+1;
        else
            k_nbr = k_cur-1;
            
        // wait until result of neighboring group is available
        cv_active.wait(ul, [k_nbr, this] {
            return this->flag_grp[k_nbr] != -1;
            });
        int up_nbr = flag_grp[k_nbr];
        
        // synchronize with neighbor
        node.wait(rank, up, up_nbr);

        // all groups have synchronized, mark tree as ready
        // todo: support more than one level (-> stack)

    }
};

int main() {
    return 0;
}
