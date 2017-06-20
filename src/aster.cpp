#ifndef STAND_ALONE
#include <Rcpp.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <queue>
//#include <unordered_map>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <climits>
#include <cmath>
#include <cfloat>

#define unordered_map map

#ifdef STAND_ALONE // pure C++ (currently not work)
#define cerr std::cerr
#else // Rcpp
#include "parent_set.h"
#define cerr Rcerr
using namespace Rcpp;
#endif

double Bayes_score(IntegerMatrix T, int m, int proc=0, double s=0, int n=0, int ss=1);
double bound(IntegerMatrix T, int m, int proc=0, int n=0, int ss=1);

//typedef long long int int64_t;
//typedef unsinged long long int uint64_t;
typedef uint64_t vset;
typedef double score_t;

const uint64_t ONE_LLU = 1;
//const score_t INF = 99999999.0;


class OrderNode {
private:
    vset vs_;
    score_t score_;
    score_t estimated_score_;
    int parent_pos_; // parent_pos_ == -1 means no parent

public:
    // for the initial node
    OrderNode(score_t estimated_score) :
        vs_(0), score_(0.0), estimated_score_(estimated_score), parent_pos_(-1) { }

    OrderNode(vset vs, score_t score, score_t estimated_score, int parent_pos) :
        vs_(vs), score_(score), estimated_score_(estimated_score), parent_pos_(parent_pos) { }

    vset getVset() const
    {
        return vs_;
    }

    score_t getScore() const
    {
        return score_;
    }

    score_t getEstimatedScore() const
    {
        return estimated_score_;
    }

    int getParentPos() const
    {
        return parent_pos_;
    }

    void setScore(score_t score)
    {
        score_ = score;
    }

    void setEstimatedScore(score_t score)
    {
        estimated_score_ = score;
    }

    void setParentPos(int parent_pos)
    {
        parent_pos_ = parent_pos;
    }

    bool isGoal(int n)
    {
        return vs_ == ((ONE_LLU << n) - 1);
    }

    std::string toString(int n)
    {
        std::ostringstream oss;
        for (int i = 0; i < n; ++i) {
            oss << (((vs_ >> i) & ONE_LLU) != 0 ? '1' : '0');
        }
        oss << ", " << score_ << ", " << estimated_score_ << ", " << parent_pos_;
        return oss.str();
    }
};

bool operator< (const OrderNode& node1, const OrderNode& node2)
{
    return node1.getEstimatedScore() < node2.getEstimatedScore();
}

bool operator> (const OrderNode& node1, const OrderNode& node2)
{
    return node1.getEstimatedScore() > node2.getEstimatedScore();
}


class ParentSet {
private:
    std::vector<std::unordered_map<vset, vset> > y_maps_;
    std::vector<std::unordered_map<vset, score_t> > z_maps_;
    
    const int sign_;
    const double s_;
    const int n_;
    const int ss_;

public:
    // If is_flip is true, the sign of scores is flipped.
    ParentSet(bool is_flip, double s=0, int n=0, int ss=0) : sign_(is_flip ? -1 : 1),
                                                             s_(s),
                                                             n_(n),
                                                             ss_(ss) { }

    void computeParentSet(NumericMatrix matrix, int tree_width = 0, int proc = 0)
    {
        int n = matrix.ncol();

        y_maps_.resize(n);
        z_maps_.resize(n);

        for (int i = 0; i < n; ++i) {
            std::vector<IntegerVector> ws;
            std::vector<IntegerVector> ys;

            DataFrame df = parent(matrix, i, tree_width, proc);

            for (int j = 0; j < n - 1; ++j) {
                char ss[32];
                sprintf(ss, "w.%d", j + 1); // 1 origin for R
                ws.push_back(df[ss]); // obtain j-th column of w
                sprintf(ss, "y.%d", j + 1); // 1 origin for R
                ys.push_back(df[ss]); // obtain j-th column of y
            }
            NumericVector z = df["z"];

            for (size_t j = 0; j < ws.size(); ++j) {
                IntegerVector vv = ws[j];
                cerr << "ws[" << j << "]: ";
                for (int k = 0; k < vv.length(); ++k) {
                    cerr << vv[k] << ", ";
                }
                cerr << "\n";
            }

            for (size_t j = 0; j < ys.size(); ++j) {
                IntegerVector vv = ys[j];
                cerr << "ys[" << j << "]: ";
                for (int k = 0; k < vv.length(); ++k) {
                    cerr << vv[k] << ", ";
                }
                cerr << "\n";
            }

            cerr << "z = ";
            for (int k = 0; k < z.length(); ++k) {
                cerr << z[k] << ", ";
            }
            cerr << "\n";

            for (int j = 0; j < z.length(); ++j) {
                vset vsw = BitvecToVset(ws, j, n, i);
                vset vsy = BitvecToVset(ys, j, n, i);
                y_maps_[i][vsw] = vsy;
                z_maps_[i][vsw] = sign_ * static_cast<score_t>(z[j]);
            }
        }
    }

    void computeParentSetEx(NumericMatrix matrix, int tree_width = 0, int proc = 0)
    {
        int ncol = matrix.ncol();

        y_maps_.resize(ncol);
        z_maps_.resize(ncol);

        std::vector<int> kind_vec;
        std::vector<int> bit_pos_vec;

        IntegerMatrix imat = normalize_table(matrix, &kind_vec, &bit_pos_vec);

        for (int i = 0; i < ncol; ++i) {
            parent_ex(imat, i, tree_width, proc,
                      &y_maps_[i], &z_maps_[i], kind_vec, bit_pos_vec);
        }
    }

    score_t getBestScore(int x, vset u_vset) const
    {
        try {
            score_t v = z_maps_.at(x).at(u_vset);
            return sign_ * v;
        } catch (...) {
            stop("not found %d at z_maps[%d]", u_vset, x);
        }
    }

    vset getBestParent(int x, vset u_vset) const
    {
        try {
            vset vs = y_maps_.at(x).at(u_vset);
            return vs;
        } catch (...) {
            stop("not found %d at y_maps[%d]", u_vset, x);
        }
    }

private:
    vset BitvecToVset(const std::vector<IntegerVector>& vec, int j, int n, int p)
    {
        vset vs = 0;
        int s = 0;
        for (int i = 0; i < n - 1; ++i) {
            if (i == p) {
                ++s;
            }
            if (vec[i][j] != 0) {
                vs |= (ONE_LLU << s);
            }
            ++s;
        }
        return vs;
    }

    int zerocount(uint64_t x)
    {
        x = x | (x >> 1);
        x = x | (x >> 2);
        x = x | (x >> 4);
        x = x | (x >> 8);
        x = x | (x >> 16);
        x = x | (x >> 32);
        return 64 - __builtin_popcountll( ~x );
    }

    bool next_permutation(uint64_t* perm, int n)
    {
        uint64_t& v = *perm;
        int p;
        for (p = 0; p < n; ++p) {
            if ((v & (ONE_LLU << (n - p - 1))) == 0) {
                break;
            }
        }
        int q = p;
        for ( ; p < n; ++p) {
            if ((v & (ONE_LLU << (n - p - 1))) != 0) {
                break;
            }
        }
        if (p == n) { // perm is 1...10...0
            return false;
        }
        v &= ~(((ONE_LLU << (p + 1)) - 1) << (n - p - 1));
        v |= (((ONE_LLU << (q + 1)) - 1) << (n - p));
        return true;
    }

    IntegerMatrix normalize_table(NumericMatrix& matrix,
                                  std::vector<int>* kind_vec,
                                  std::vector<int>* bit_pos_vec)
    {
        int nrow = matrix.nrow();
        int ncol = matrix.ncol();
        IntegerMatrix ret(nrow, ncol);

        kind_vec->resize(ncol);
        bit_pos_vec->resize(ncol + 1);
        (*bit_pos_vec)[0] = 0;

        for (int j = 0; j < ncol; ++j) {
            std::map<int, int> val_map;
            int c = 0;
            for (int i = 0; i < nrow; ++i) {
                if (val_map.find(matrix(i, j)) == val_map.end()) { // not found
                    val_map[matrix(i, j)] = c;
                    ++c;
                }
                ret(i, j) = val_map[matrix(i, j)];
            }
            (*kind_vec)[j] = c;
            (*bit_pos_vec)[j + 1] = (*bit_pos_vec)[j] + zerocount(c - 1);
            if ((*bit_pos_vec)[j + 1] > 64) {
                stop("S cannot exceed 64!");
            }
        }
        return ret;
    }

    DataFrame parent_ex(IntegerMatrix matrix, int h, int tw, int proc,
                   std::map<uint64_t, uint64_t>* y_map,
                   std::map<uint64_t, double>* z_map,
                   const std::vector<int>& kind_vec,
                   const std::vector<int>& bit_pos_vec)
    {
        int col = matrix.ncol();
        if (tw == 0) {
            tw = col - 1;
        }
        std::map<uint64_t, uint64_t>& y = *y_map;
        std::map<uint64_t, double>& z = *z_map;

        std::map<uint64_t, bool> x;

        uint64_t mask = (ONE_LLU << h) - 1;

        for (int i = 0; i <= tw; ++i) {
            uint64_t perm = (ONE_LLU << i) - 1;
            do {
                uint64_t ex_perm = (perm & mask) | ((perm & ~mask) << 1);
                x[ex_perm] = false;
                z[ex_perm] = -100000000.0;
                for (int j = 0; j < col; ++j) {
                    if ((ex_perm & (ONE_LLU << j)) != 0) {
                        uint64_t child = (ex_perm & ~(ONE_LLU << j));
                        if (x[child]) {
                            x[ex_perm] = true;
                        }
                        if (z[child] > z[ex_perm]) {
                            y[ex_perm] = y[child];
                            z[ex_perm] = z[child];
                        }
                    }
                }
                if (!x[ex_perm]) {
                    int m = kind_vec[h];
                    IntegerMatrix T = fftable_ex(matrix, m, h, ex_perm, bit_pos_vec);
                    if (z[ex_perm] > bound(T, m, proc, n_, ss_)) {
                        x[ex_perm] = true;
                    } else {
                        double s = Bayes_score(T, m, proc, s_, n_, ss_);
                        if (s > z[ex_perm]) {
                            y[ex_perm] = ex_perm;
                            z[ex_perm] = s;
                        }
                    }
                }
            } while (next_permutation(&perm, col - 1));
        }
        return DataFrame();
    }

    IntegerMatrix fftable_ex(IntegerMatrix& matrix, int w, int h,
                             uint64_t children, const std::vector<int>& bit_pos_vec)
    {
        int nrow = matrix.nrow();
        int ncol = matrix.ncol();
        IntegerMatrix T(nrow, w);
        std::map<uint64_t, std::map<int, int> > counter;

        for (int i = 0; i < nrow; ++i) {
            uint64_t b = 0;
            for (int j = 0; j < ncol; ++j) {
                if ((children & (ONE_LLU << j)) != 0) {
                    b |= (static_cast<uint64_t>(matrix(i, j)) << bit_pos_vec[j]);
                }
            }
            int val = matrix(i, h);
            counter[b][val] += 1;
        }
        int row = 0;
        std::map<uint64_t, std::map<int, int> >::iterator itor1 = counter.begin();
        while (itor1 != counter.end()) {
            int col = 0;
            std::map<int, int>::iterator itor2 = itor1->second.begin();
            while (itor2 != itor1->second.end()) {
                T(row, col) = itor2->second;
                ++col;
                ++itor2;
            }
            ++row;
            ++itor1;
        }
        return T(Range(0, row - 1), _);
    }
};


class ASterQueue {
private:
    const int n_;
    const int tree_width_;
    std::vector<vset> heap_; // for the priority queue
    std::unordered_map<vset, OrderNode*> node_map_;
    std::unordered_map<vset, int> pos_map_;
    std::vector<std::unordered_map<vset, score_t> > cache_map_;
    const ParentSet& parent_set_;

public:
    ASterQueue(int n, int tree_width, const ParentSet& parent_set) : n_(n),
        tree_width_((tree_width != 0) ? tree_width : n),
        cache_map_(n),
        parent_set_(parent_set)
    { }

    void addNode(OrderNode* node)
    {
        node_map_[node->getVset()] = node;
        heap_.push_back(node->getVset());
        int index = static_cast<int>(heap_.size()) - 1;
        pos_map_[node->getVset()] = index;
        upHeap(index);
    }

    void addOrUpdateNode(vset vs, score_t score, score_t estimated_score, int parent_pos)
    {
        std::unordered_map<vset, OrderNode*>::iterator itor = node_map_.find(vs);
        if (itor != node_map_.end()) { // found
            if (estimated_score < itor->second->getEstimatedScore()) { // update
                itor->second->setScore(score);
                itor->second->setEstimatedScore(estimated_score);
                itor->second->setParentPos(parent_pos);
                upHeap(pos_map_[vs]);
            }
        } else { // not found
            OrderNode* node = new OrderNode(vs, score, estimated_score, parent_pos);
            node_map_[node->getVset()] = node;
            heap_.push_back(node->getVset());
            int index = static_cast<int>(heap_.size()) - 1;
            pos_map_[node->getVset()] = index;
            upHeap(index);
        }
    }

    OrderNode* pop()
    {
        assert(heap_.size() > 0);
        int last_index = static_cast<int>(heap_.size()) - 1;
        swapPos(0, last_index);

        vset vs = heap_.back();
        OrderNode* node = node_map_[vs];
        heap_.pop_back();
        pos_map_.erase(vs); // vs of the keys in node_map_ is not erased
                            // because it will be used later
        downHeap(0);
        return node;
    }

    size_t size() const
    {
        return heap_.size();
    }

    bool empty() const
    {
        return heap_.empty();
    }

    OrderNode* getNode(vset vs) const
    {
        return node_map_.at(vs);
    }

    score_t getBestScore(int x, vset u_vset)
    {
        if (__builtin_popcountll(u_vset) <= tree_width_) {
            score_t v = parent_set_.getBestScore(x, u_vset);
            return v;
        } else {
            score_t cand_score = -99999999.9;
            for (int y = 0; y < n_; ++x) {
                if (((u_vset >> y) & ONE_LLU) != 0) {
                    vset u_minus_y_vset = u_vset & ~(ONE_LLU << y);
                    score_t v;
                    std::unordered_map<vset, score_t>::iterator itor =
                        cache_map_[x].find(u_minus_y_vset);
                    if (itor != cache_map_[x].end()) { // element found
                        v = itor->second;
                    } else {
                        v = getBestScore(x, u_minus_y_vset);
                    }
                    if (v > cand_score) {
                        cand_score = v;
                    }
                }
            }
            cache_map_[x][u_vset] = cand_score;
            return cand_score;
        }
    }

    void print()
    {
        std::vector<vset>::iterator itor = heap_.begin();
        while (itor != heap_.end()) {
            cerr << *itor << ", ";
            ++itor;
        }
        cerr << "\n";
    }

    void printheap(int n)
    {
        std::vector<vset>::iterator itor = heap_.begin();
        while (itor != heap_.end()) {
            cerr << *itor << ": " << node_map_[*itor]->toString(n) << std::endl;
            ++itor;
        }
    }

    void test()
    {
        addNode(new OrderNode(0, 0, 5.0, 0));
        addNode(new OrderNode(1, 0, 2.0, 0));
        addNode(new OrderNode(2, 0, 3.0, 0));
        addNode(new OrderNode(3, 0, 6.0, 0));
        addNode(new OrderNode(4, 0, 1.0, 0));
        addNode(new OrderNode(5, 0, 4.0, 0));
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        addNode(new OrderNode(6, 0, 3.0, 0));
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        printheap(5);cerr << "---" << std::endl;
        OrderNode* node = new OrderNode(7, 0, 2.0, 0);
        addNode(node);
        addNode(new OrderNode(8, 0, 8.0, 0));
        addNode(new OrderNode(9, 0, 7.0, 0));
        printheap(5);cerr << "---" << std::endl;
        node->setEstimatedScore(10.0);
        update(node->getVset());
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        printheap(5);
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
        cerr << pop()->getEstimatedScore() << ", " << std::endl;
    }

private:
    void swapPos(int index1, int index2)
    {
        pos_map_[heap_[index1]] = index2;
        pos_map_[heap_[index2]] = index1;
        std::swap(heap_[index1], heap_[index2]);
    }

    void upHeap(int index)
    {
        if (index > 0) {
            int parent = (index - 1) / 2;
            if (node_map_[heap_[parent]]->getEstimatedScore() >
                node_map_[heap_[index]]->getEstimatedScore()) {
                swapPos(parent, index);
                upHeap(parent);
            }
        }
    }

    void downHeap(int index)
    {
        int child = -1;
        int n = static_cast<int>(heap_.size());

        if (index * 2 + 2 < n) {
            child = (node_map_[heap_[index * 2 + 1]]->getEstimatedScore() <
                     node_map_[heap_[index * 2 + 2]]->getEstimatedScore() ?
                     index * 2 + 1 : index * 2 + 2);
        } else if (index * 2 + 2 == n) {
            child = index * 2 + 1;
        }
        if (child >= 0) { // a child exists
            if (node_map_[heap_[child]]->getEstimatedScore() <
                node_map_[heap_[index]]->getEstimatedScore()) {
                swapPos(child, index);
                downHeap(child);
            }
        }
    }

    void update(vset vs)
    {
        int index = pos_map_[vs];
        int parent = (index - 1) / 2;
        if (node_map_[heap_[parent]]->getEstimatedScore() >
            node_map_[heap_[index]]->getEstimatedScore()) {
            upHeap(index);
        } else {
            downHeap(index);
        }
    }
};


std::vector<int> extractResult(int n, const ASterQueue& queue, const ParentSet& parent_set)
{
    std::vector<int> result;
    vset node = (ONE_LLU << n) - 1;
    while (node > 0) {
        int x = queue.getNode(node)->getParentPos();
        vset u = (ONE_LLU << x);
        node &= ~u;
        vset c = parent_set.getBestParent(x, node);
        for (int i = 0; c != 0; ++i) {
            //if (i != x) { // skip if i == x
                if ((c & ONE_LLU) != 0) { // make arc "i -> x"
                    result.push_back(i);
                    result.push_back(x);
                }
                c >>= 1;
            //}
        }
    }
    return result;
}

score_t heuristics(vset vs)
{
    return 0;
}

// Todo: add const property to ParentSet
std::vector<int> runAster(NumericMatrix& matrix, const ParentSet& parent_set, int tree_width)
{
    const int n = matrix.ncol();
    ASterQueue queue(n, tree_width, parent_set);

    score_t initial_estimated_score = heuristics(0);
    queue.addNode(new OrderNode(initial_estimated_score)); // initial node

    while (!queue.empty()) {

        // pop the node having the minimum estimated score
        // u must not be deleted because it will be used later.
        OrderNode* u = queue.pop();

        if (u->isGoal(n)) { // reach the goal
            cerr << "goal: " << u->toString(n) << "\n";
            return extractResult(n, queue, parent_set);
        }
        vset u_vset = u->getVset();
        for (int x = 0; x < n; ++x) {
            if (((u_vset >> x) & ONE_LLU) == 0) { // for each i in V minus U
                vset x_vset = (ONE_LLU << x);
                vset ux = (u_vset | x_vset);
                score_t g_new = queue.getBestScore(x, u_vset) + u->getScore();
                queue.addOrUpdateNode(ux, g_new, g_new + heuristics(ux), x);
            }
        }
    }
    return std::vector<int>();
}


// [[Rcpp::export]]
NumericVector aster_cpp(NumericMatrix matrix, int tree_width = 0, int proc=1, double s=0, int n=0, int ss=0)
{
    if (matrix.ncol() >= 64) {
        stop("The data with more than 64 variables is not supported.");
    }

    // give true to the argument because of flipping the sign
    // for solving the maximization problem.
    ParentSet parent_set(true, s, n, ss);
    //parent_set.computeParentSet(matrix, tree_width, proc);
    parent_set.computeParentSetEx(matrix, tree_width, proc);

    // double start_time = sec();
    // cerr << "Algorithm Start\n";
    std::vector<int> result = runAster(matrix, parent_set, tree_width);
    // double end_time = sec();
    // cerr << "time = " << end_time - start_time << " sec." << std::endl;

    NumericVector v;
    for (size_t i = 0; i < result.size(); ++i) {
        v.push_back(result[i]);
    }
    return v;
}

