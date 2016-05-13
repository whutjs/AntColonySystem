#ifndef ANTS_H
#define ANTS_H

#include <vector>
#include <set>
#include <map>

/* add a small constant to avoid division by zero if a distance is 
zero */
#define EPSILON            0.00000000000001
#define MAX_NEIGHBOR	   8		// Ã¿¸öµã×î¶àÖ»ÓÐ¸ö³ö¶È£¬Ò²¾ÍÊÇ8¸öÁÚ¾Ó
#define MAX_NODES 615		  // max no. of nodes
#define INFTY 100000
#define MAX_REQ_NODES 55		  // ×î¶àµÄÖÐ¼ä½Úµã
#define MAX_EDGES 4810			// ×î¶àµÄ±ßÊý


// ±£´æ½ÚµãµÄÐÅÏ¢
struct NodeInfo {
	unsigned short indegree;		// Èë¶È
	unsigned short outdegree;		// ³ö¶È
	std::vector<unsigned short> indegree_list;		// Èë¶È±í
	std::vector<unsigned short> outdegree_list;		// ³ö¶È±í
	unsigned short out_include_reqnode_cnt;			// ³ö¶È°üº¬±Ø¾­µãµÄ¸öÊý
	NodeInfo() {
		indegree = outdegree = out_include_reqnode_cnt = 0;
	}
};

struct Edge{
    unsigned short linkid;			// the max. of linkid is 4800
    unsigned short cost;			
    Edge() {
        cost = 10000;
    }
};

struct Ant{
   std::vector<unsigned short>  tour;	// the path
   std::vector<unsigned short> req_node_idx;
   unsigned short curNode;			// the node this ant currently in
   unsigned short requireCnt;		// how many require nodes this ant has passed
   int  tour_length;
   bool visited[MAX_NODES];
   bool stop;				// whether this ant is in stop state
   bool isForward;			// ÊÇ·ñÊÇ´ÓÆðµãµ½ÖÕµãµÄËÑË÷£¬false±íÊ¾´ÓÖÕµãµ½ÆðµãµÄËÑË÷
};
/* for local search */
// 保存的是将非对称转换为对称之后的节点信息，在local_search里面会用到
extern std::map<unsigned short, NodeInfo> node_map_for_ls;
// 保存的是影子节点的偏移，=nodes.size()+1
extern int SHADOW_OFFSET;
/* 保存的是虚拟节点的偏移（也就是到该节点权重都为0的点），=SHADOW_OFFSET+nodes.size()+1
	虚拟节点的node_id = 某非必经点node_id+VIRT_OFFSET（1200）
	而虚拟节点的影子节点的id=虚拟节点的node_id+SHADOW_OFFSET（600）
	所以node_id的最大值是:599(正常节点id),1199(正常节点的影子节点的最大id),1799(虚拟节点的最大id)
	2399(虚拟节点的影子节点的最大id)，所以WeightMatrix必须是600*4
	*/
extern int VIRT_OFFSET;
struct WeightMatrix{
	int cost;
	// 使用类的原因是为了使用构造函数
	WeightMatrix(){ cost = 10000;}
};
extern WeightMatrix weightMatrix[MAX_NODES*4][MAX_NODES*4];
// 保存虚拟节点
extern std::vector<unsigned short> virt_nodes;

/* end local search */

// key=nodeid,  value=info
extern std::map<unsigned short, NodeInfo> nodeInfoMap;

// ±£´æÁË:Æðµãµ½ËùÓÐ±Ø¾­µãµÄÂ·¾¶£¬±Ø¾­µãÁ½Á½Ö®¼ä×î¶ÌÂ·¾¶£¬±Ø¾­µãµ½ÖÕµã×î¶ÌÂ·¾¶
extern std::vector<std::vector<unsigned short> > total_shortest_path;

// key=Ò»¸ö·Ç±Ø¾­µãv£¬value=Ò»¸ö±Ø¾­µãv', ¶øv'Ö»ÓÐvÄÜµ½´ï£¬ËùÒÔÔÚvÑ¡ÔñÏÂÒ»¸ö½ÚµãµÄÊ±ºò±ØÐëÒªµ½´ïv'
extern std::map<unsigned short, unsigned short> mustGoToReqNode;
extern std::set<unsigned short> notFoundReqNode;		// ±£´æ»¹Ã»ÕÒµ½µÄÖÐ¼ä½Úµã
extern std::set<unsigned short> notFoundReqNodeIndegreeList;		// »¹Ã»ÕÒµ½µÄÖÐ¼ä½ÚµãµÄÈë¶È

extern int reqNodeMaxOutDegree, nodeMaxIndegree;		// ±Ø¾­µãµÄ×î´ó³ö¶ÈºÍËùÓÐ½áµãµÄ×î´óÈë¶È
//extern double best_ant_weight;			// ×î¼ÑÂ·¾¶ÂìÒÏ¸üÐÂµÄÈ¨ÖØ
//extern double *best_ant_weight_arrays;		// ×î¼ÑÂ·¾¶ÂìÒÏ¸üÐÂµÄÈ¡ÖµÈ¨ÖØ
//extern int best_ant_weight_arrays_length;

extern unsigned short startNode, endNode;
extern Ant *ant;      /* this (array of) struct will hold the colony */
extern Ant *best_so_far_ant;   /* struct that contains the best-so-far ant */
extern Ant *best_reverse_ant;		//	best ant for reverse direction
extern Ant* iteration_best_ant;		// ´Ë´Îµü´ú×îÓÅµÄÂìÒÏ
extern Ant *sub_opt_so_far_ant;		//	对于次最优的蚂蚁的路径，也进行奖励(差于best_so_far_ant)


extern std::map<unsigned short,double> smell[MAX_NODES];

extern double   pheromone[MAX_NODES][MAX_NODES]; /* pheromone matrix, one entry for each arc */
extern double   total[MAX_NODES][MAX_NODES];     /* combination of pheromone times heuristic information */
extern bool ignoreWeightMode;
extern Edge matrix[MAX_NODES][MAX_NODES];  // matrix[i][j] contains edge for node i->j
extern std::set<unsigned short> nodes;         // contains all nodes id
extern bool isNodesRequire[MAX_NODES];		// an array of bool indicate whether the nodeid is required(size=MAX_NODES)
extern unsigned short numOfReqNodes;		// the number of require nodes;
extern unsigned short requireNodes[MAX_REQ_NODES];

extern double  *prob_of_selection;	

extern  unsigned short n_ants;      /* number of ants */

extern double rho;           /* parameter for evaporation */
extern double alpha;         /* importance of trail */
extern double beta;          /* importance of heuristic evaluate */
extern double gama;			// for smell
extern double q_0;           /* probability of best choice in tour construction */

extern double  trail_0;         /* initial pheromone trail level in ACS  and BWAS */
extern double trail_max;
extern double trail_min;

double heuristic(const unsigned short &from, const unsigned short &to);
/* Pheromone manipulation etc. */

void init_pheromone_trails (const double &initial_trail);

void global_update_pheromone (Ant &k);

void compute_total_information( );

/* Ants' solution construction */

void ant_empty_memory(Ant &k);

// replace the ant in a startNode
void re_place_ant(Ant &k);
void place_ant_to(Ant &k, unsigned short node);


bool neighbour_choose_and_move_to_next( Ant &k);;
/* Auxiliary procedures related to ants */

int find_best ();
int find_most_midnode();
bool copy_from_to(Ant *a1, Ant *a2);

void allocate_ants();


void global_acs_pheromone_update(Ant *k );

void local_acs_pheromone_update(Ant *k );

int compute_tour_length(std::vector<unsigned short> &path);

int compute_reverse_tour_length(std::vector<unsigned short> &path);


#endif
