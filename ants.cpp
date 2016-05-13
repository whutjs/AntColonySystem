#include "stdafx.h"
#include "ants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "direct_search.h"

using namespace std;

/* for local search */
// 保存的是将非对称转换为对称之后的节点信息，在local_search里面会用到
map<unsigned short, NodeInfo> node_map_for_ls;
// 保存的是影子节点的偏移，=nodes.size()+1
int SHADOW_OFFSET = 600;
// 保存的是虚拟节点的偏移（也就是到该节点权重都为0的点），=SHADOW_OFFSET+nodes.size()+1
int VIRT_OFFSET = 1200;
// 保存的是ATSP->TSP之后的权重信息
WeightMatrix weightMatrix[MAX_NODES*4][MAX_NODES*4];
vector<unsigned short> virt_nodes;
/* end local search */

unsigned short startNode, endNode;
Ant *ant;      /* this (array of) struct will hold the colony */
Ant *best_so_far_ant;   /* struct that contains the best-so-far ant */
Ant *best_reverse_ant;
Ant* iteration_best_ant;		// ´Ë´Îµü´ú×îÓÅµÄÂìÒÏ
Ant *sub_opt_so_far_ant;		//	对于次最优的蚂蚁的路径，也进行奖励(差于best_so_far_ant)

map<unsigned short, NodeInfo> nodeInfoMap;


// ±£´æÁË:Æðµãµ½ËùÓÐ±Ø¾­µãµÄÂ·¾¶£¬±Ø¾­µãÁ½Á½Ö®¼ä×î¶ÌÂ·¾¶£¬±Ø¾­µãµ½ÖÕµã×î¶ÌÂ·¾¶
vector<std::vector<unsigned short> > total_shortest_path;

map<unsigned short, unsigned short> mustGoToReqNode;
int reqNodeMaxOutDegree, nodeMaxIndegree;	// ±Ø¾­µãµÄ×î´ó³ö¶ÈºÍËùÓÐ½áµãµÄ×î´óÈë¶È
set<unsigned short> notFoundReqNode;
set<unsigned short> notFoundReqNodeIndegreeList;		// »¹Ã»ÕÒµ½µÄÖÐ¼ä½ÚµãµÄÈë¶È

double  pheromone[MAX_NODES][MAX_NODES]; /* pheromone matrix, one entry for each arc */
double  total[MAX_NODES][MAX_NODES];     /* combination of pheromone times heuristic information */
/* ±£´æÒ»¸öµãµÄÆøÎ¶,Ò»¹²×î¶àÓÐ600¸öµã£¬Ã¿¸öµã±£´æµÄÊÇÄ³¸ö±Ø¾­µã¶Ô¸ÃµãµÄÆøÎ¶
	µÄÓ°Ïì£¬Èçmap = smell[0]±íÊ¾µÄÊÇ½Úµã0µÄÆøÎ¶ÐÅÏ¢£¬ÆäÖÐmap[node_id]±íÊ¾µÄÊÇ
	±Ø¾­µãnode_id¶Ô¸ÃµãµÄÆøÎ¶µÄÓ°Ïì 
*/
map<unsigned short,double> smell[MAX_NODES];				

Edge matrix[MAX_NODES][MAX_NODES];  // 只保存正常节点的边信息
vector<unsigned short> bestPath;

set<unsigned short> nodes;         // contains all nodes id
bool isNodesRequire[MAX_NODES];		// an array of bool indicate whether the nodeid is required(size=MAX_NODES)
unsigned short requireNodes[MAX_REQ_NODES];

unsigned short numOfReqNodes;		// the number of require nodes;
double  *prob_of_selection;	
bool ignoreWeightMode;
unsigned short n_ants;      /* number of ants */
double rho;           /* parameter for evaporation */
double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
double gama;
double q_0;           /* probability of best choice in tour construction */
double  trail_0;         /* initial pheromone trail level in ACS */
double trail_max;
double trail_min;


double heuristic(const unsigned short &from, const unsigned short &to){
	if(ignoreWeightMode) {
		return 1;		
	}else{
		return 1.0 / ((double) matrix[from][to].cost);
	}
}
void init_pheromone_trails(const double &initial_trail ) {
    //printf(" init trails with %.5f\n",initial_trail);
    /* Initialize pheromone trails */
	for(set<unsigned short>::iterator from = nodes.begin(); from != nodes.end(); from++) {
		for(set<unsigned short>::iterator to = nodes.begin(); to != nodes.end(); to++) {
			pheromone[*from][*to] = initial_trail;			
			total[*from][*to] = initial_trail;
		}
	}
	
	/*
	if (ignoreWeightMode == false) {   
        // ¶ÔÒªÇóµÄ½ÚµãºÍÖÕµã¼ÓÉÏÒ»¶¨µÄÈ¨ÖØ
	for(set<unsigned short>::iterator from = nodes.begin(); from != nodes.end(); from++) {
		bool flag = false;
		if(isNodesRequire[*from]) {
			flag = true;
		}
		for(int i = 0; i < adjList[*from].size(); i++) {
			int to = adjList[*from][i];
			if(to == endNode || isNodesRequire[to]){
				//pheromone[*from][to] += 1.0/matrix[*from][to].cost;
					pheromone[*from][to] += initial_trail;
				if(flag) {
					pheromone[*from][to] += initial_trail;
				}
			}			
		}
	}
	}
	*/
	
	
}

// reinforces edges used in ant a's solution
void global_update_pheromone(Ant &k){     
   // printf("global pheromone update\n");
    double d_tau = 1.0 / (double) k.tour_length;
	unsigned short from = 0, to = 0;
	for(int i = 0; i < (k.tour.size() - 1); i++) {
		from  = k.tour[i];
		to = k.tour[i+1];
		pheromone[from][to] += d_tau;
	}   
}

// calculates heuristic info times pheromone for each arc
void compute_total_information( ) {     
    //printf("compute total information\n");	
	for (set<unsigned short>::iterator from = nodes.begin(); from != nodes.end(); from++) {
		for (set<unsigned short>::iterator to = nodes.begin(); to != nodes.end(); to++){
			if(*from == *to || matrix[*from][*to].cost == 10000) continue;
			total[*from][*to] = pow(pheromone[*from][*to], alpha) * 
								pow(heuristic(*from,*to),beta);				
		}
    }
}

void ant_empty_memory(Ant &k) {	
	for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		k.visited[*it] = false;
	}
	k.stop = false;
	k.requireCnt = 0;
	k.tour.clear();
	k.tour_length = 0;
}

void re_place_ant(Ant &k) {
	k.curNode = startNode;
	k.tour.clear();
	k.tour.push_back(startNode);
	k.visited[startNode] = true;
	k.isForward = true;
	// in case the require nodes include the startNode
	if(isNodesRequire[startNode]) {
		k.requireCnt = 1;
	}
}

void place_ant_to(Ant &k, unsigned short node){
	k.curNode = node;
	k.tour.clear();
	k.tour.push_back(node);
	k.visited[node] = true;
	// in case the require nodes include the startNode
	if(isNodesRequire[node]) {
		k.requireCnt++;
	}
}
/*
	return: true if this ant stops moving, otherwise false.
*/
bool neighbour_choose_and_move_to_next(Ant &k) {  
    double   rnd, partial_sum = 0., sum_prob = 0.0;
	rand();	
	bool is_forward = k.isForward;
	unsigned short currentNode = k.curNode, nextNode = 0; 	
	if (is_forward) {
		// Ö»ÓÐÇ°ÏòËÑË÷²ÅÐèÒª£¬ÒòÎª´ÓÖÕµãµ½ÆðµãµÄËÑË÷ÊÇÓÃÈë¶ÈËÑË÷£¬Èç¹ûÖ»ÓÐÒ»¸öµã¾ÍÖ»»áÈ¥ÄÇ
		map<unsigned short, unsigned short>::iterator mustIt = mustGoToReqNode.find(currentNode);
		if (mustIt != mustGoToReqNode.end()) {
			//printf("mustGoToReqNode\n");
			bool stop = false;
			// ´Óµ±Ç°½Úµã±ØÐëÒªµ½µÄ½Úµã
			nextNode = mustIt->second;
			if (k.visited[nextNode] == false) {
				//printf("from %d must go to %d\n",currentNode, nextNode);
				k.tour.push_back(nextNode);
				k.visited[nextNode] = true;
				if (isNodesRequire[nextNode]) {
					k.requireCnt++;		
					k.req_node_idx.push_back(k.tour.size()-1);
				}
				if (nextNode == endNode || nextNode == startNode) {
					k.stop = true;
					stop = true;
				}

				k.curNode = nextNode;
				return stop;
			}
		}
	}
	if (is_forward) {
		// µ±Ö»Ê£1¸öÒÔÏÂÃ»ÓÐ¾­¹ýµÄ½ÚµãÊ±£¬Ò»Óöµ½¾Íµ½
		if (notFoundReqNode.size() == 1) {
		//	printf("notFoundReqNode.size() == 1\n");
			if (notFoundReqNodeIndegreeList.find(currentNode) != notFoundReqNodeIndegreeList.end()) {
				double p = 1.0 / notFoundReqNodeIndegreeList.size();
				double x = ((double)rand() / RAND_MAX);
				// ÒÔÒ»¶¨µÄ¸ÅÂÊÑ¡Ôñ±ØÐëµ½´ïµ±Ç°½Úµã
				if (x <= p) {
					bool stop = false;
					// ´Óµ±Ç°½Úµã±ØÐëÒªµ½Õâ¸ö½Úµã
					nextNode = *(notFoundReqNode.begin());
					//printf("from %d must go to %d,p=%f x=%f\n", currentNode, nextNode, p, x);
					if (k.visited[nextNode] == false) {
						//printf("from %d must go to %d\n",currentNode, nextNode);
						k.tour.push_back(nextNode);
						k.visited[nextNode] = true;
						if (isNodesRequire[nextNode]) {
							k.requireCnt++;							
						}
						if (nextNode == endNode || nextNode == startNode) {
							k.stop = true;
							stop = true;
						}
						k.curNode = nextNode;
						return stop;
					}
				}
			}
		}
	}
	const NodeInfo& nodeInfo = nodeInfoMap[k.curNode];
	
	// ÓÃÀ´±éÀúÁÚ½Ó±íµÄµü´úÆ÷
	vector<unsigned short>::const_iterator adjlist_begin;
	vector<unsigned short>::const_iterator adjlist_end;
	int nn_ants = 0;
	if(is_forward) {
		// Ç°ÏòËÑË÷ÓÃµÄÊÇ³ö¶È
		nn_ants = nodeInfo.outdegree;
		adjlist_begin = nodeInfo.outdegree_list.begin();
		adjlist_end = nodeInfo.outdegree_list.end();
	}else{
		nn_ants = nodeInfo.indegree;
		adjlist_begin = nodeInfo.indegree_list.begin();
		adjlist_end = nodeInfo.indegree_list.end();
	}
	 double  *prob_ptr = prob_of_selection;
	//double  prob_ptr[10];
	prob_ptr[nn_ants] = HUGE_VAL;
	
	int idx = 0;
	bool hasEndNode = false;
	for(; adjlist_begin != adjlist_end; adjlist_begin++) 
	{
		unsigned short to = *adjlist_begin;
		if(k.visited[to] || currentNode == to) {
			prob_ptr[idx++] = 0;
		}		
		else if((to == endNode || to == startNode) && k.requireCnt < numOfReqNodes){
			// ÖÕµãµ¥¶À´¦Àí
			// Èç¹û»¹Ã»·ÃÎÊ¹»×ã¹»µÄµã£¬ÔòÏÈ²»ÈÃ·ÃÎÊÖÕµã
			hasEndNode = true;
			prob_ptr[idx++] = 0;
		}		
		else{
			// TODO£º¶ÔÃ¿Ö»ÂìÒÏÀ´Ëµ£¬¼ÓÈëµÄÆøÎ¶ÊÇ²»Ò»ÑùµÄ
			double smell_weight = 0.0001;
			//if(ignoreWeightMode && nodes.size() != 598) {				
			if(nodes.size() != 598) {				
				//unsigned short link_id = matrix[currentNode][to].linkid;	
				//printf("cur_node=%d to=%d \n", currentNode, to);
				const map<unsigned short, double>& smell_info = smell[to];
				map<unsigned short, double>::const_iterator m_it = smell_info.begin();
				for(; m_it != smell_info.end(); m_it++) {
					unsigned short smell_req_node = m_it->first;
					if(isNodesRequire[smell_req_node] && k.visited[smell_req_node] == false) {
						// find the larger one
						if(m_it->second > smell_weight)
						    smell_weight = m_it->second;
					}
				}
			}
			prob_ptr[idx] = total[currentNode][to];
			//if(ignoreWeightMode == false && nodes.size() != 598) {				
			if(nodes.size() != 598) {				
				prob_ptr[idx] *= pow(smell_weight, gama);	
			//	prob_ptr[idx] *= smell_weight;
			}
			sum_prob += prob_ptr[idx];
			idx++;
		}
	}
    if (sum_prob <= 0.0) {
		
		if(hasEndNode && k.requireCnt < numOfReqNodes){
			// ³ýÁËÖÕµãÒÔÍâÃ»ÓÐ±ðµÄµã¿ÉÒÔ·ÃÎÊÁË
			unsigned short tnode = startNode;
			if(is_forward) {
				tnode = endNode;
			}
			k.tour.push_back(tnode);
			k.visited[tnode] = true;
			k.stop =  true;
			k.curNode = tnode;
			return true;
		}else{
			// All cities from the adjList are visited 
			//printf("sum_prob <= 0.0 ant in curNode=%d stop\n", k.curNode);
			k.stop = true;
			return true;
		}
		
		//k.stop = true;
		//return true;
    }     
    else { 		
		bool stop = false;
		int to = 0;
		rnd = ((double) rand() / RAND_MAX);
		rnd *= sum_prob;
		to = 0;
		partial_sum = prob_ptr[to];	
		while (partial_sum <= rnd) {
			to++;
			if(to >= nn_ants){
				to = nn_ants-1;
				break;
			}
			partial_sum += prob_ptr[to];			
		}
		if(is_forward) {
			nextNode = nodeInfo.outdegree_list[to];
		}else{
			nextNode = nodeInfo.indegree_list[to];
		}
		k.tour.push_back(nextNode);
		k.visited[nextNode] = true;
		if(isNodesRequire[nextNode]) {
			k.requireCnt++;
		}		
		if(nextNode == endNode || nextNode == startNode) {
			//printf("nextNode == endNode ant in curNode=%d stop\n", k.curNode);
			k.stop = true;
			stop = true;
		}
		k.curNode = nextNode;
		return stop;
    }
}

int find_most_midnode() {
	int   k, k_min;
	int maxNodeCnt = 0;
	k_min = -1;
	double density = 0;
	int length = INFTY;
	for (k = 0; k < n_ants; k++) {
		const Ant& ak = ant[k];
		if (ak.requireCnt >= numOfReqNodes && (ak.curNode == startNode || ak.curNode == endNode)
			&& ak.tour_length < length) {
			k_min = k;
			length = ak.tour_length;
			continue;
		}
		double d = ak.requireCnt / (double)ak.tour.size();
		if (ak.requireCnt > maxNodeCnt && d > density) {
			maxNodeCnt = ak.requireCnt;
			k_min = k;
		}
	}
	if(k_min != -1) {
		if (ant[k_min].requireCnt > iteration_best_ant->requireCnt) {
			copy_from_to(&ant[k_min], iteration_best_ant);				
		}
	}
	return k_min;
}

/* ÕÒµ½ËùÓÐÂìÒÏÖÐ£¬Â·¾¶×î¶ÌµÄ */
 int find_best( ) {
     int   min = INFTY;
     int   k, k_min = -1;
	int maxreqcnt_reverse = 0, k_min_reverse = -1;
    for( k = 0 ; k < n_ants ; k++ ) {		
		if((ant[k].curNode == endNode || ant[k].curNode == startNode) && ant[k].requireCnt >= numOfReqNodes 
			&& ant[k].tour_length < min ) {
			min = ant[k].tour_length;
			k_min = k;
		}
    }
	if(k_min == -1) {
		// ÕÒµ½¾­¹ý×î¶à±Ø¾­µãµÄ·´ÏòÂìÒÏ
		 for( k = 0 ; k < n_ants ; k++ ) {		
			 if(ant[k].isForward == false) {
				 if(ant[k].requireCnt > maxreqcnt_reverse) {
					 maxreqcnt_reverse = ant[k].requireCnt;
					 k_min_reverse = k;
				 }
			 }
		 }
		 if(k_min_reverse != -1) {
			 copy_from_to(best_reverse_ant, &ant[k_min_reverse]);
		 }
	}else{
		// 找到一个最优解后，再找一个次优解,此时min=best_so_far_ant.tourLength
		int sub_opt_min = INFTY, sub_opt_k = -1;
		 for( k = 0 ; k < n_ants ; k++ ) {		
			 const Ant& ak = ant[k];
			 unsigned short req_cnt = ak.requireCnt;
			 if(req_cnt < numOfReqNodes) continue;
			 if((ak.curNode != endNode && ak.curNode != startNode)) continue;
			 int tour_length = ak.tour_length;			
		      if(tour_length > min && tour_length < sub_opt_min) {
				sub_opt_min = ant[k].tour_length;
				sub_opt_k = k;
			}
		}
		 if(sub_opt_k != -1) {
			// printf("found a sub_opt_min, length=%d\n", sub_opt_min);
			// for(int i = 0; i < ant[sub_opt_k].tour.size(); i++) printf("%d ", ant[sub_opt_k].tour[i]);
			// printf("\n");
			 if(sub_opt_min < sub_opt_so_far_ant->tour_length){
				 copy_from_to(&ant[sub_opt_k], sub_opt_so_far_ant);
			 }
		 }
	}
    return k_min;
}

 /* 
	copy solution from ant a1 into ant a2
	return: false if failed
 */
 bool copy_from_to(Ant *a1, Ant *a2) {
	 if(a1 == a2 || a1 == NULL || a2 == NULL) return false;
	 a2->req_node_idx.clear();
	 a2->req_node_idx.insert(a2->req_node_idx.end(), 
		 a1->req_node_idx.begin(), a1->req_node_idx.end());
	 a2->curNode = a1->curNode;
	 a2->requireCnt = a1->requireCnt;
	 a2->stop = a1->stop;
	 a2->tour_length = a1->tour_length;
	 a2->tour.clear();
	 a2->tour.insert(a2->tour.begin(), a1->tour.begin(), a1->tour.end());
	 a2->isForward = a1->isForward;
		for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++)
		{
			a2->visited[*it] = a1->visited[*it];
		}
	 return true;
	 /*
	 for(int i = 0; i < MAX_NODES; i++) {
		 a2.visited[i] = a1.visited[i];
	 }
	 */
 }

 void allocate_ants() {	     
    if((ant = new Ant[n_ants]) == NULL){
		printf("Out of memory, exit.");
		exit(1);
    }
	for(int i = 0; i < n_ants; i++) {
		ant[i].curNode = startNode;
		ant[i].requireCnt = 0;
		ant[i].stop = false;
		ant[i].tour_length = 0;
		for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++) {
			ant[i].visited[*it] = false;
        }
	}
	if((best_so_far_ant = new Ant) == NULL){
		printf("Out of memory, exit.");
		exit(1);
    }
	 best_so_far_ant->tour_length = INFTY;
	 best_so_far_ant->requireCnt = 0;
	 best_so_far_ant->curNode = 900;
	 best_so_far_ant->isForward = true;
		 
	if((sub_opt_so_far_ant = new Ant) == NULL){
		printf("Out of memory, exit.");
		exit(1);
    }
	 sub_opt_so_far_ant->tour_length = INFTY;
	 sub_opt_so_far_ant->requireCnt = 0;
	 sub_opt_so_far_ant->curNode = 900;
	 sub_opt_so_far_ant->isForward = true;

	if((best_reverse_ant = new Ant) == NULL){
		printf("Out of memory, exit.");
		exit(1);
    }
	 best_reverse_ant->tour_length = INFTY;
	 best_reverse_ant->requireCnt = 0;
	 best_reverse_ant->curNode = 900;
	 best_reverse_ant->isForward = false;
	 
    if((iteration_best_ant = new Ant) == NULL){
		printf("Out of memory, exit.");
		exit(1);
     }
	 iteration_best_ant->tour_length = INFTY;
	 iteration_best_ant->requireCnt = 0;
	 iteration_best_ant->curNode = 900;
	 iteration_best_ant->isForward = false;

	 int max_node_degree = reqNodeMaxOutDegree;
	 if(nodeMaxIndegree > max_node_degree) {
		 max_node_degree = nodeMaxIndegree;
	 }
	if ((prob_of_selection = new double[max_node_degree + 2]) == NULL) {
		printf("Out of memory, exit.");
		exit(1);
    }   
     prob_of_selection[max_node_degree] = HUGE_VAL;
 }

 void global_acs_pheromone_update(Ant *k) {  
	 if(k->tour_length == INFTY) return;
	// printf("acs_global_update, best length:%d\n", best_so_far_ant->tour_length);
	double weight = n_ants;
	//double weight = 1.0;
//	printf("best_ant_weight=%f\n", best_ant_weight);
    double d_tau = weight / (double) k->tour_length;
	unsigned short from, to;
	if(k->isForward) {
		for (int i = 0 ; i < k->tour.size()-1; i++ ) {
			from = k->tour[i];
			to = k->tour[i+1];
			if(pheromone[from][to] <= trail_max){
				pheromone[from][to] = (1. - rho) * pheromone[from][to] + rho * d_tau;	
				// ½±Àø¼ÆËãÌ«¶à£¬·ÅÔÚÒ»ÆðÖØÐÂ¼ÆËãtotal
				// total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);		
			}
		}
	}else{
		for (int i = k->tour.size()-1 ; i > 0 ; i--) {
			from = k->tour[i];
			to = k->tour[i-1];
			if(pheromone[from][to] <= trail_max){
				pheromone[from][to] = (1. - rho) * pheromone[from][to] + rho * d_tau;	
				// ½±Àø¼ÆËãÌ«¶à£¬·ÅÔÚÒ»ÆðÖØÐÂ¼ÆËãtotal
				// total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);		
			}else{
				//pheromone[from][to] = trail_max;
			}
		}
	}
}

 /* 
	removes some pheromone on edge just passed by the ant
	Note:should do experiments to find the best value of parameter xi. 
 */
 void local_acs_pheromone_update(Ant *k) {  
	 if(k == NULL) return;
	 if(k->tour.size() < 2){
		 printf("in local_acs_pheromone_update, error!tour size is less than 2!\n");
		 return;
	 }
	 int tourSize = k->tour.size();
	 unsigned short from = 0, to = 0;	
	 if(k->isForward) {
		 from = k->tour[tourSize-2];
		 to = k->tour[tourSize-1];
	 }else{
		 from = k->tour[tourSize-1];
		 to = k->tour[tourSize-2];
	 }
	/* still additional parameter has to be introduced */
	 double xi = 0.1;
	pheromone[from][to] = (1. - xi) * pheromone[from][to] + xi * trail_0;    
	total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);
}

 int compute_tour_length(vector<unsigned short> &path){
	int size = path.size();
    int tour_length = 0;  
	unsigned short from = 0, to = 0;
    for (int i = 0 ; i < size-1 ; i++ ) {
		from = path[i];
		to = path[i+1];
		tour_length += matrix[from][to].cost;
    }
    return tour_length;
}

int compute_reverse_tour_length(vector<unsigned short> &path){
	int size = path.size();
    int tour_length = 0;  
	unsigned short from = 0, to = 0;
    for (int i = size-1 ; i > 0 ; i--) {
		from = path[i];
		to = path[i-1];
		tour_length += matrix[from][to].cost;
    }
    return tour_length;
}

