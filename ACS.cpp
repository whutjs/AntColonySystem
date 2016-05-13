#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <set>
#include <math.h>
#include <assert.h>
#include "ACS.h"
#include "time.h"
#include "stdlib.h"
#include "direct_search.h"
#include "lib_record.h"
#include <string.h>

using namespace std;


int found_best;		// ֒սخԅޢքּպˇ֚һՎּպ
int n_tours;		// ٹլנʙٶ·޶èһһּ֨ˇࠉѐ·޶é
int iteration;		// ּպ݆˽Ƿ
int max_tours = 0;		// һՎtryՊѭٹլנʙٶtours
double max_time;		// Պѭִѐքʱݤ
int max_tries;
double punish_fac_for_zero;		// הԚޭڽ0ٶҪȳ֣ք·޶քԍףӲؓ
double punish_fac_for_notend;		// הԚûԐսկָ֨ו֣քԍףӲؓ

double alphas_weight[10] = {1, 1, 2, 1, 3, 2, 2, 1, 1, 2};		// ԃ4̬֯ټтхϢ̘alpha
double betas_weight[10] =  {2, 3, 2, 4, 3, 3, 4, 4, 5, 5};		// ԃ4̬֯ټтхϢ̘beta

double alphas[10] = { 2, 2, 2, 2, 2, 2, 2, 2, 1, 1 };             // 用来动态更新信息素alpha
double betas[10] = { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 };             // 用来动态更新信息素beta

//double rhos[10] = {0.1, 0.2, 0.3, 0.4, 0.9, 0.6, 0.7, 0.8, 0.5, 0.65};
double rhos[10] = {0.1, 0.2, 0.3, 0.4, 0.9, 0.6, 0.7, 0.8, 0.5, 0.65};
int rhos_length = 10;
int rhoIdx = 9;

int lastBest = 0;			// ʏһՎiterationքخݑֵìɧڻlѸ5Վּպּˇһҹìղהrho޸ѐҤۻ
int sameBestCnt = 0;		// ɧڻlѸ5Վּպخ׌·޶ּˇһҹìղהrho޸ѐҤۻ

int checkNodeNotFound();
void sort_nodemap_for_ls();
void add_to_node_map_for_ls(const unsigned short& from, const unsigned short& to, 
							const unsigned short& cost);
void print_nodemap_for_ls();
// 在local_search里面用到的visited数组
bool *visited_for_ls;

bool termination_condition( ) {
	double passTime = elapsed_time( VIRTUAL );
	return  ((n_tours >= max_tours) || sameBestCnt >= ((max_tours>>1)+1) || ( passTime >= max_time));  
}

void construct_solutions(){
     int k;        /* counter variable */    
   // printf("construct solutions for all ants\n");

    /* Mark all cities as unvisited */
    for ( k = 0 ; k < n_ants ; k++) {
		ant_empty_memory(ant[k]);
    }        
    /* Place the ants on same initial city */
   // for ( k = 0 ; k < n_ants ; k++ )
	//	re_place_ant(ant[k]);
	// 使用两个方向的蚂蚁
	for ( k = 0 ; k < n_ants ; k++ ){
		if(k % 2 == 0) {
			place_ant_to(ant[k], startNode);
			ant[k].isForward = true;
		}else{
			place_ant_to(ant[k], endNode);
			ant[k].isForward = false;
		}
	}
	int stopAnts = 0; 
    while ( stopAnts < n_ants ) {		
		for ( k = 0 ; k < n_ants ; k++ ) {
			if(ant[k].stop) continue;
			bool stop = neighbour_choose_and_move_to_next(ant[k]);			
			//printf("ant[%d] choose:%d\n", k, ant[k].curNode);
			//local_acs_pheromone_update( &ant[k]);
			if(stop) {
				stopAnts++;
			}
		}		
    }  
    for ( k = 0 ; k < n_ants ; k++ ) {
		if(ant[k].isForward) {
			ant[k].tour_length = compute_tour_length(ant[k].tour);
		}else{
			ant[k].tour_length = compute_reverse_tour_length(ant[k].tour);		
		}
    }
    n_tours += 1;
}

// 在找到一条可行解之后，用dijkstra找到的最短路指导寻路
void rewardForDijkstraShortPath() {
	if(total_shortest_path.size() == 0) return;
	double weight = 1.0;	
    double d_tau = weight / (double) best_so_far_ant->tour_length;
	vector<vector<unsigned short> >::const_iterator pathSetIt = total_shortest_path.begin();
	vector<vector<unsigned short> >::const_iterator pathSetEnd = total_shortest_path.end();
	while(pathSetIt != pathSetEnd) {
		const vector<unsigned short>& path = *pathSetIt;
		if(path.size() <= 1){
			pathSetIt++;
			continue;
		}
		// 路径是从后往前的
		for(int i = path.size()-1; i > 0; i--) {
			unsigned short from = path[i], to = path[i-1];
			// TODO:加多少是个问题。。。
			pheromone[from][to] = (1. - rho) * pheromone[from][to] + rho*d_tau;
			if(pheromone[from][to] >= trail_max) {
				pheromone[from][to] = trail_max;
			}
		}
		pathSetIt++;
	}
}

// 强化最短路径
void rewardForShortPath() {
	if(nodes.size() == 598) return;
	if(total_shortest_path.size() == 0) return;
	// 对初始最短路径加上一定的权重
	vector<vector<unsigned short> >::const_iterator pathSetIt = total_shortest_path.begin();
	vector<vector<unsigned short> >::const_iterator pathSetEnd = total_shortest_path.end();
	while(pathSetIt != pathSetEnd) {
		const vector<unsigned short>& path = *pathSetIt;
		if(path.size() <= 1){
			pathSetIt++;
			continue;
		}
		// 路径是从后往前的
		for(int i = path.size()-1; i > 0; i--) {
			unsigned short from = path[i], to = path[i-1];
			pheromone[from][to] += trail_min * (1./rho);
			//pheromone[from][to] += trail_min;
			if(pheromone[from][to] >= trail_max) {
				pheromone[from][to] = trail_max;
			}
		}
		pathSetIt++;
	}
}
// 奖励此次迭代最优的蚂蚁
void rewardForIterationBest(){
	if (iteration_best_ant->tour_length == INFTY) return;
	const Ant* k = iteration_best_ant;
	//printf("acs specific: global pheromone update\n");	
	double weight = k->requireCnt;	
	double d_tau = 0;
	if (ignoreWeightMode == false) {
		d_tau = weight / (double)k->tour_length;
	}
	else {
		d_tau = weight / (double)k->tour.size() * (weight / numOfReqNodes);
	}
	unsigned short from, to;
	if (k->isForward) {
		for (int i = 0; i < k->tour.size() - 1; i++) {
			from = k->tour[i];
			to = k->tour[i + 1];
			pheromone[from][to] = (1. - rho) * pheromone[from][to] + rho * d_tau;
			if(pheromone[from][to] >= trail_max) {
				pheromone[from][to] = trail_max;
			}
			// 奖励计算太多，放在一起重新计算total
			//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from, to), beta);
		}
	}
	else {
		for (int i = k->tour.size() - 1; i > 0; i--) {
			from = k->tour[i];
			to = k->tour[i - 1];
			pheromone[from][to] = (1. - rho) * pheromone[from][to] + rho * d_tau;
			if(pheromone[from][to] >= trail_max) {
				pheromone[from][to] = trail_max;
			}
			// 奖励计算太多，放在一起重新计算total
			//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from, to), beta);
		}
	}
}

void punishAndReward() {
	double bestR = best_so_far_ant->requireCnt / (double)best_so_far_ant->tour.size();
	for(int i = 0; i < n_ants; i++) {
		if(ant[i].isForward == false)continue;
		if(ant[i].requireCnt == 0) {
			//printf("punish ant:%d for zero\n", i);
			unsigned short from = 0, to = 0;
			for(int j = 0; j < ant[i].tour.size()-1; j++) {
				from = ant[i].tour[j];
				to = ant[i].tour[j+1];
				if(pheromone[from][to] > trail_min) {
						pheromone[from][to] *= punish_fac_for_zero;
						//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);
				}else{
					pheromone[from][to] = trail_min;
				}									
			}
		}
		if(ant[i].curNode != endNode) {
			//printf("punish ant:%d for not end node\n", i);
			unsigned short from = 0, to = 0;
			for(int j = 0; j < ant[i].tour.size()-1; j++) {
				from = ant[i].tour[j];
				to = ant[i].tour[j+1];
				if(isNodesRequire[to] == false) {
					if(pheromone[from][to] > trail_min) {
						pheromone[from][to] *= punish_fac_for_notend;
						//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);
					}else{
						pheromone[from][to] = trail_min;
					}					
				}
			}
		}
		if(best_so_far_ant->requireCnt < numOfReqNodes && ant[i].curNode == endNode && ant[i].requireCnt > 0) {
			//printf("reward for the ant:%d right path\n", i);
			unsigned short from = 0, to = 0;
			//double weight = ant[i].requireCnt;
			double weight = 1.0;
			int antTourLength = ant[i].tour_length;
			if(antTourLength < best_so_far_ant->tour_length) {
				antTourLength = best_so_far_ant->tour_length;
			}
			for(int j = 0; j < ant[i].tour.size()-1; j++) {
				from = ant[i].tour[j];
				to = ant[i].tour[j+1];
				if(to == endNode || isNodesRequire[to]) {
					if(pheromone[from][to] < trail_max){
						pheromone[from][to] += (weight / antTourLength) * rho;
						//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);
					}else{
					//	pheromone[from][to] = trail_max;
					}
				}
			}
		}
	}
	// reward for best reverse ant
	if(best_so_far_ant->requireCnt < numOfReqNodes) {
			//printf("reward for the ant:%d right path\n", i);
			const Ant& r_ant = *best_reverse_ant;
			unsigned short from = 0, to = 0;			
			double weight = 1.0;
			int antTourLength = r_ant.tour_length;
			if(antTourLength < best_so_far_ant->tour_length) {
				antTourLength = best_so_far_ant->tour_length;
			}
			int size = r_ant.tour.size();
			for(int j = size - 1; j > 0; j--) {
				from = r_ant.tour[j];
				to = r_ant.tour[j-1];
				if(pheromone[from][to] < trail_max){
					pheromone[from][to] += (weight / antTourLength) * rho;
					//total[from][to] = pow(pheromone[from][to], alpha) * pow(heuristic(from,to),beta);
				}else{
					//pheromone[from][to] = trail_max;
				}
			}
		}
	// 奖励计算太多，放在一起重新计算total
	// compute_total_information();
}

void update_statistics() {
	if( best_so_far_ant->requireCnt >= numOfReqNodes) {
		double p_x = exp(log(0.05)/nodes.size()); 
		trail_min = 1. * (1. - p_x) / (p_x * (double)((8) / 2));
		//trail_min = 0.01;
		trail_max = 1. / (best_so_far_ant->tour_length );
		trail_0 = trail_max;
		trail_min = trail_max * trail_min; 
	}
	int iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */
	find_most_midnode();
	if(iteration_best_ant == -1) {
		//printf("solution not found\n");
		return;
	}	
    if ( ant[iteration_best_ant].tour_length < best_so_far_ant->tour_length ) {	
		copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
		found_best = iteration;		
    }
	//printf("found best:%d\n", ant[iteration_best_ant].tour_length);
	//for(int i = 0; i < best_so_far_ant->tour.size(); i++) printf("%d ", best_so_far_ant->tour[i]);
	//printf("\n");
}

void acs_global_update() {	
	//printf("acs_global_update for length:%d\n", best_so_far_ant->tour_length);
    global_acs_pheromone_update( best_so_far_ant );
	//printf("acs_global_update for subopt length:%d\n", sub_opt_so_far_ant->tour_length);
    //global_acs_pheromone_update( sub_opt_so_far_ant );
}
void pheromone_trail_update() {   
    //acs_global_update();
	// 信息素挥发
	set<unsigned short>::const_iterator fromIt = nodes.begin();
	set<unsigned short>::const_iterator toIt = nodes.begin();
	for( ; fromIt != nodes.end(); fromIt++) {
		for( ; toIt != nodes.end(); toIt++) {
			unsigned short from = *fromIt, to = *toIt;
			if(from == to) continue;
			if(matrix[from][to].cost <= 20) {
				 pheromone[from][to] = (1 - rho) * pheromone[from][to];
				 if(pheromone[from][to] <= trail_min) {
					 pheromone[from][to] = trail_min;
				 }
			}
		}
	}
	//printf("pheromone_trail_update\n");
	acs_global_update();
	// 奖励计算太多，放在一起重新计算total
	//compute_total_information();
}


/* 
   occasionally compute some statistics and check whether algorithm 
   is converged
*/
void search_control_and_statistics() {
	if(iteration % 5 == 0) {
               // printf("iteration=%d\n", iteration);   
		rhoIdx--;
		if(rhoIdx <= 0) rhoIdx = rhos_length-1;
		rho = rhos[rhoIdx];
		//printf("new rho=%f\n", rho);
	}
	
	int bestLength = best_so_far_ant->tour_length;
	if(bestLength != INFTY && bestLength == lastBest) {
      //  printf("bestLength=%d lastBest=%d\n", bestLength, lastBest);
		sameBestCnt++;
	}else{
		sameBestCnt = 0;
	}
	if(bestLength < INFTY) {
		lastBest = bestLength;
	}

}

/* initilialize variables appropriately when starting a trial */
void init_try(const int& try_cnt){
    //best_ant_weight = best_ant_weight_arrays[ (try_cnt%best_ant_weight_arrays_length)];
	if(ignoreWeightMode){
		alpha = alphas[try_cnt];
		beta = betas[try_cnt];
	}else{
		alpha = alphas_weight[try_cnt];
		beta = betas_weight[try_cnt];
	}
    //printf("init try,alpha=%.5f beta=%.5f\n", alpha, beta);  
    rhoIdx = rhos_length - 1;
    rho = rhos[rhoIdx];
    sameBestCnt = 0;
    /* Initialize variables concerning statistics etc. */
    n_tours      = 0;
    iteration    = 1;       
    found_best   = 0;
    init_pheromone_trails(trail_0);
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();
    //printf("trail_max=%f, trail_min=%f\n", trail_max, trail_min);
}

void set_default_parameters() {    
    ignoreWeightMode = true;
    n_ants         = nodes.size()/1.5 + 0.5;    /* number of ants (-1 means instance size) */    
	// n_ants         = 30;
    alpha          = 2.0;
    beta           = 2.0;
	gama = 1.0;
    rhoIdx         =  rhos_length - 1;
    rho            = rhos[rhoIdx];
    q_0            = 0.1;
  //  best_ant_weight_arrays_length = 3;
   // best_ant_weight_arrays = new double[best_ant_weight_arrays_length];
   // best_ant_weight_arrays[0] = 1.0;
  //  best_ant_weight_arrays[1] = n_ants;
   // best_ant_weight_arrays[2] = numOfReqNodes;
    if(nodes.size() <= 5) {
		max_tries = 1;
        max_tours = 10;
    }else if(nodes.size() <= 10){
		  max_tries = 2;
          max_tours = 30;
    }else{
		max_tries = 9;
        max_tours = 90;
    }    
    //seed           = ( int) time(NULL);
    max_time       = 9.92;
   // optimal        = 1;
	//trail_0 = 1.0 / INFTY;
	trail_max = 1.0 / MAX_NODES;
	trail_min = trail_max / ( 2. * nodes.size());
	trail_0 = trail_max;

	punish_fac_for_zero = 0.4;
	punish_fac_for_notend = 0.7;

	//printf("ants number=%d\n", n_ants);	
}

/* 将原图转换为atsp图（非必经点全连通)
	也就是新增nodes.size()-2-numOfReqNodes+1个节点（除去必经点，起点和终点，
	剩余的就是非必经点，再加1个点）
*/
void trans_graph_to_atsp() {
	SHADOW_OFFSET = 600;
	VIRT_OFFSET = 1200;
	/* 虚拟节点按endNode->第一个非必经点->第二个非必经点..的顺序
	 把原图串联起来。因此在virt_nodes里面保存的也是这个顺序的virtnode
	 */
	// 首先加的虚拟节点A=endNode+VIRT_OFFSET,用来把终点连接起来
	unsigned short cost_to_vnode = 0;
	unsigned short vnode1_id = endNode + VIRT_OFFSET;
	add_to_node_map_for_ls(endNode, vnode1_id, cost_to_vnode);
	// 保存的是所有非必经点+起点
	vector<unsigned short> node_to_link;
	for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++) {
		unsigned short nreq_node = *it;
		if(!isNodesRequire[nreq_node] && nreq_node != startNode && nreq_node != endNode ) {
			node_to_link.push_back(nreq_node);
		}
	}
	// 最后加入起点
	node_to_link.push_back(startNode);
	int node_to_link_idx = 0;
	int node_to_link_size = node_to_link.size();
	queue<unsigned short> vnode_queue;
	vnode_queue.push(vnode1_id);
	/* 接下来，将每个在队列头里的虚拟节点vnode_n出队 
	 然后将vnode_n连接到下一个还没连接的非必经点v。同时生成一个
	 新的虚拟节点vnode_n+1，将vnode_n连接到vnode_n+1,同时将v也连接到vnode_n+1
	*/
	while(node_to_link_idx < node_to_link_size) {
		unsigned short cur_node_to_link = node_to_link[node_to_link_idx];
		unsigned short cur_vnode = vnode_queue.front();
		virt_nodes.push_back(cur_vnode);
		vnode_queue.pop();
		// link cur_vnode to next_node_to_link
		add_to_node_map_for_ls(cur_vnode, cur_node_to_link, cost_to_vnode);
		if(cur_node_to_link == startNode) {
			// 已经link到终点了，break
			break;
		}
		// 然后生成一个新的vnode
		unsigned short next_vnode = cur_node_to_link + VIRT_OFFSET;
		// 将cur_vnode->next_vnode, cur_node_to_link->next_vnode
		add_to_node_map_for_ls(cur_vnode, next_vnode, cost_to_vnode);
		add_to_node_map_for_ls(cur_node_to_link, next_vnode, cost_to_vnode);
		// 入队
		vnode_queue.push(next_vnode);
		node_to_link_idx++;
	}

	// 为了防止某些点在node_map_for_ls少加了出度点，这里对所有的点都加上自己作为出度
	for(set<unsigned short>::iterator nit = nodes.begin(); nit != nodes.end();nit++) {
		unsigned short A = *nit;
		unsigned short A_p = 600 + A;
		map<unsigned short, NodeInfo>::iterator AIt = node_map_for_ls.find(A);
		if (AIt != node_map_for_ls.end())
		{
			// 将A_p加到A的出度表里
			// 出度+1
			AIt->second.outdegree += 1;
			AIt->second.outdegree_list.push_back(A_p);		
		}
		else {
			//printf("from not found, outdegree = 1\n");
			NodeInfo nodeInfo;
			nodeInfo.outdegree = 1;
			nodeInfo.indegree = 0;			
			nodeInfo.outdegree_list.push_back(A_p);
			node_map_for_ls[A] = nodeInfo;
		}

		map<unsigned short, NodeInfo>::iterator Ap_It = node_map_for_ls.find(A_p);
		if (Ap_It != node_map_for_ls.end())
		{
			// 将A加到A_p的出度表里
			// 出度+1
			Ap_It->second.outdegree += 1;
			Ap_It->second.outdegree_list.push_back(A);		
		}
		else {
			//printf("from not found, outdegree = 1\n");
			NodeInfo nodeInfo;
			nodeInfo.outdegree = 1;
			nodeInfo.indegree = 0;			
			nodeInfo.outdegree_list.push_back(A);
			node_map_for_ls[A_p] = nodeInfo;
		}
	}

	// 为了防止虚拟点在node_map_for_ls少加了自己的影子节点作为出度点，这里对所有的点都加上自己作为出度
	for(int i = 0; i < virt_nodes.size(); i++) {		
		unsigned short A = virt_nodes[i];
		unsigned short A_p = 600 + A;
		map<unsigned short, NodeInfo>::iterator AIt = node_map_for_ls.find(A);
		if (AIt != node_map_for_ls.end())
		{
			// 将A_p加到A的出度表里
			// 出度+1
			AIt->second.outdegree += 1;
			AIt->second.outdegree_list.push_back(A_p);		
		}
		else {
			//printf("from not found, outdegree = 1\n");
			NodeInfo nodeInfo;
			nodeInfo.outdegree = 1;
			nodeInfo.indegree = 0;			
			nodeInfo.outdegree_list.push_back(A_p);
			node_map_for_ls[A] = nodeInfo;
		}

		map<unsigned short, NodeInfo>::iterator Ap_It = node_map_for_ls.find(A_p);
		if (Ap_It != node_map_for_ls.end())
		{
			// 将A加到A_p的出度表里
			// 出度+1
			Ap_It->second.outdegree += 1;
			Ap_It->second.outdegree_list.push_back(A);		
		}
		else {
			//printf("from not found, outdegree = 1\n");
			NodeInfo nodeInfo;
			nodeInfo.outdegree = 1;
			nodeInfo.indegree = 0;			
			nodeInfo.outdegree_list.push_back(A);
			node_map_for_ls[A_p] = nodeInfo;
		}
	}
	sort_nodemap_for_ls();
	print_nodemap_for_ls();
}

void initProgram() {			
	visited_for_ls = new bool[MAX_NODES*4];
	set_default_parameters();
	allocate_ants();	
	shortest_path_between_reqnode(total_shortest_path);	
	printf("total_shortest_path size=%d\n", total_shortest_path.size());
	init_smell();
	// init local search
	/*
	for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++){
			for(set<unsigned short>::iterator to = nodes.begin(); to != nodes.end(); to++){
				if(*to == *it) continue;
				if(matrix[*it][*to].cost == 0) continue;

				printf("%d -> %d smell=%f\n", *it, *to, smell[matrix[*it][*to].linkid]);
			}
		}
	*/
}

bool construct_first_tour() {
	ignoreWeightMode = true;
	bool found = false;
	int loopTime = 0;
	double local_max_time = 9.8;
	while (best_so_far_ant != NULL && best_so_far_ant->requireCnt < numOfReqNodes
		&& elapsed_time(VIRTUAL) < local_max_time) {			
		for (int n_try = 0; found == false && n_try < max_tries; n_try++) {		
			if(n_try > 0) {
				// 增加迭代次数
				max_tours += 10;
				 if(max_tours >= 300) max_tours = 300;
			}
			init_try(n_try);
			while (n_tours < max_tours && elapsed_time(VIRTUAL) < local_max_time) {
				/*
				printf("pheromenon:\n");
				for(set<unsigned short>::const_iterator fromIt = nodes.begin(); fromIt != nodes.end(); fromIt++){
					for(set<unsigned short>::const_iterator toIt = nodes.begin(); toIt != nodes.end(); toIt++){
						if(matrix[*fromIt][*toIt].cost > 0)
							printf("(%d,%d)=%f ", *fromIt, *toIt, pheromone[*fromIt][*toIt]);
					}
					printf("\n");
				}
				*/
				rewardForShortPath();
				construct_solutions();
				/*
				printf("construct_solutions:----------------\n\n");
				for(int i = 0; i < n_ants; i++) {
					printf("ant[%d], length=%d stop:%d requireCnt=%d path:\n", i,  ant[i].tour_length, (ant[i].stop==true), ant[i].requireCnt);
					for(int j = 0; j < ant[i].tour.size(); j++) {
						printf("%d ", ant[i].tour[j]);
					}
					printf("\n\n");
				}
				*/
				update_statistics();
				if (best_so_far_ant->tour_length < INFTY) {
					found = true;
					break;
				}		
				pheromone_trail_update();
				punishAndReward();
				rewardForIterationBest();
				search_control_and_statistics();
				iteration++;
				checkNodeNotFound();
				// 奖励计算太多，放在一起重新计算total
				compute_total_information();
			}
		}
		loopTime++;
	}
	return found;
}


// 检查还有多少个中间节点没有找到
int checkNodeNotFound() {
	int cnt = 0;
	notFoundReqNode.clear();
	notFoundReqNodeIndegreeList.clear();
	if (best_so_far_ant->requireCnt < numOfReqNodes) {
		map<unsigned short, bool> foundNodes;
		for (int i = 0; i < numOfReqNodes; i++) {
			foundNodes[requireNodes[i]] = false;
		}
		for (int i = 0; i < best_so_far_ant->tour.size(); i++) {
			unsigned short node = best_so_far_ant->tour[i];
			if (isNodesRequire[node]) {
				foundNodes[node] = true;
			}
		}
		map<unsigned short, bool>::iterator it;
		for (it = foundNodes.begin(); it != foundNodes.end(); ++it) {
			if (it->second == false) {
				cnt++;
				notFoundReqNode.insert(it->first);
				const NodeInfo& nodeInfo = nodeInfoMap[it->first];
				notFoundReqNodeIndegreeList.insert(nodeInfo.indegree_list.begin(), nodeInfo.indegree_list.end());
			}
		}
	}
	return cnt;
}

/* is_forward:true表示路径是正向 */
void local_search(vector<unsigned short>& tour,  vector<int>& tsp_tour,
				  bool is_forward) {
	if(tour.size() == 0) return;
	for(int i = 0; i < MAX_NODES*4; i++)
		visited_for_ls[i] = false;
	// 对tour改造成对称TSP的形式
	
	int idx = 0, end_idx = tour.size();
	int step = 1;
	if(is_forward == false) {
		idx = tour.size() - 1;
		end_idx = -1;
		step = -1;
	}
	while(idx != end_idx) {
		unsigned short A = tour[idx];
		unsigned short A_p = tour[idx] + 600;
		tsp_tour.push_back(A);
		tsp_tour.push_back(A_p);
		idx += step;
	}
	for(int i = 0; i < tsp_tour.size(); i++) {
		unsigned short node = tsp_tour[i];
		visited_for_ls[node] = true;
	}
	// 然后将剩下的没有被访问过的非必经点加入路径
	unsigned short end_node = tsp_tour[tsp_tour.size()-1];
	queue<unsigned short> Q;
	Q.push(end_node);
	//visited_for_ls[end_node] = false;
	while(Q.empty() == false) {
		int cur_node_id = Q.front();
		Q.pop();
		if(cur_node_id != end_node)
			tsp_tour.push_back(cur_node_id);
		if(cur_node_id == startNode) break;
		visited_for_ls[cur_node_id] = true;
		/* 从当前节点寻找下一个可能的节点,对于剩下的节点，每个节点
		最多只可能有两个出度:(1)对于虚拟节点，一个出度点是非必经点，
		另一个出度点是相邻的虚拟节点；(2)对于非必经点，只有一个出度点：虚拟节点
		*/
		const NodeInfo& cur_node = node_map_for_ls[cur_node_id];
		bool is_vnode = (cur_node_id >= VIRT_OFFSET ? true : false);
		if(!is_vnode) {
			/* 不是虚拟节点，因为出度都已经按降序排过序，所以选择的是第一个还没被访问过的点。 
				理由:假如是A,那么A的出度有A'(-10000)，V(虚拟节点，0），其他正常节点(>0)。
				那么如果A的影子节点A'没有被访问，那就应该到A'。如果已经访问了，那就应该
				是虚拟节点。
			*/
			int next_node = -1;
			for(int i = 0; i < cur_node.outdegree_list.size(); i++) {
				int  tnode = cur_node.outdegree_list[i];
				if(visited_for_ls[tnode] == false){
					next_node = tnode;
					break;				
				}
			}
			if(next_node != -1)
				Q.push(next_node);
		}else{
			// 是虚拟节点，有两个选择：(1)下一个虚拟节点,(2)非必经点
			// 虚拟节点的出度中，也是有:到自己影子节点A'是-10000，到非必经点0
			// 但如果有非必经点而且没被访问，那必须是该非必经点
			int next_node = -1;
			for(int i = 0; i < cur_node.outdegree_list.size(); i++) {				
				// 小于1200的是非必经点或者非必经点的影子节点
				unsigned short  tnode = cur_node.outdegree_list[i];
				if(tnode < VIRT_OFFSET && !visited_for_ls[tnode]) {
					next_node = cur_node.outdegree_list[i];
					break;				
				}
			}
			// 虚拟节点肯定会有非必经点出度，否则出错
		//	if(next_node < 0 || next_node >= 2400){
			//	printf("next_node=%d\n", next_node);
			//}
			//assert(next_node >= 0 && next_node < 2400);
			if(next_node == -1 || visited_for_ls[next_node]) {
				// 已经访问过了，那么到达下一个非访问虚拟节点
				for(int i = 0; i < cur_node.outdegree_list.size(); i++) {
					int tnode = cur_node.outdegree_list[i];
					if(tnode >= VIRT_OFFSET && !visited_for_ls[tnode]){
						next_node = tnode;
						break;				
					}
				}				
			}
			if(next_node != -1)
				Q.push(next_node);
		}
	}
	// 添加起点
	tsp_tour.push_back(tsp_tour[0]);
	return;
	int n_size = tsp_tour.size() - 1;
	//printf("before 3opt,tour_size=%d\n", n_size);
	//for(int i = 0; i < tsp_tour.size(); i++) printf("%d ", tsp_tour[i]);
	//printf("\n");	
	printf("before 3opt:length=%d\n", best_so_far_ant->tour_length);
	//printf("after 3opt:\n");
	//for(int i = 0; i < tsp_tour.size(); i++) printf("%d ", tsp_tour[i]);
	//printf("\n");
	// 再把路径转换回正确的路径
	vector<unsigned short> result;
	for(int i = 0; i < tsp_tour.size()-1; i++) {
		unsigned short node = tsp_tour[i];		
		if(node < SHADOW_OFFSET) {
			result.push_back(node);
		}
		if(node == endNode) break;
	}
	int len = compute_tour_length(result);
	printf("after 3opt:length=%d\n", len);
	//for(int i = 0; i < result.size(); i++) printf("%d ", result[i]);
	//	printf("\n");
	if(len < best_so_far_ant->tour_length) {				
		best_so_far_ant->isForward = true;		
		best_so_far_ant->tour_length = len;
		best_so_far_ant->tour.clear();
		best_so_far_ant->tour.insert(best_so_far_ant->tour.end(), result.begin(), result.end());
	}
}

// cmp_outdegree比较函数里面，比较的是哪个点到出度a,b的距离
unsigned short cmp_source_node;
/*  将node_map_for_ls里的出度表按照权重从小到大排序，方便local_search */
bool cmp_outdegree(const unsigned short &a, const unsigned short &b)
{
	int cost_a = weightMatrix[cmp_source_node][a].cost;
	int cost_b = weightMatrix[cmp_source_node][b].cost;
	if (cost_a < cost_b)
        return true;
    else return false;
}

/* 将节点加进node_map_for_ls里 */
void add_to_node_map_for_ls(const unsigned short& from, const unsigned short& to, 
							const unsigned short& cost)
{
	int A = from, B = to;
	int A_p = A + 600, B_p = B + 600;
	weightMatrix[A_p][B].cost = weightMatrix[B][A_p].cost = cost;
	// A到A'的距离是0 dist[A_p][A] = 负无穷;
	weightMatrix[A_p][A].cost = weightMatrix[A][A_p].cost = -10000;
	weightMatrix[B_p][B].cost = weightMatrix[B][B_p].cost = -10000;
	// 起点不要保存入度
	map<unsigned short, NodeInfo>::iterator BIt = node_map_for_ls.find(B);
	if (BIt != node_map_for_ls.end())
	{
		// 将A_p加到B的出度表里
		// 出度+1
		BIt->second.outdegree += 1;
		BIt->second.outdegree_list.push_back(A_p);		
	}
	else {
		//printf("from not found, outdegree = 1\n");
		NodeInfo nodeInfo;
		nodeInfo.outdegree = 1;
		nodeInfo.indegree = 0;			
		nodeInfo.outdegree_list.push_back(A_p);
		node_map_for_ls[B] = nodeInfo;
	}
	map<unsigned short, NodeInfo>::iterator Ap_It = node_map_for_ls.find(A_p);
	if (Ap_It != node_map_for_ls.end())
	{
		// 将B加到A_p的出度表里
		// 出度+1
		Ap_It->second.outdegree += 1;
		Ap_It->second.outdegree_list.push_back(B);		
	}
	else {
		//printf("from not found, outdegree = 1\n");
		NodeInfo nodeInfo;
		nodeInfo.outdegree = 1;
		nodeInfo.indegree = 0;			
		nodeInfo.outdegree_list.push_back(B);
		node_map_for_ls[A_p] = nodeInfo;
	}
}

void sort_nodemap_for_ls() {
	// 将node_map_for_ls里的所有出度排序
	for(map<unsigned short, NodeInfo>::iterator it = node_map_for_ls.begin();
		it != node_map_for_ls.end(); it++)
	{
		NodeInfo& node = (*it).second;
		cmp_source_node = (*it).first;
		sort(node.outdegree_list.begin(), node.outdegree_list.end(), cmp_outdegree);
	}
}

void print_nodemap_for_ls() {
	map<unsigned short, NodeInfo>::iterator lit = node_map_for_ls.begin();
	for(; lit != node_map_for_ls.end(); lit++) {
		unsigned short cur_node_id = lit->first;
		const NodeInfo& node = lit->second;
		printf("cur node=%d, outdegree:\n", cur_node_id);
		vector<unsigned short>::const_iterator out_it = node.outdegree_list.begin();
		for(; out_it != node.outdegree_list.end(); out_it++) {
			unsigned short out_nodeid = *out_it;
			printf("%d , cost=%d\n", out_nodeid, weightMatrix[cur_node_id][out_nodeid]);
		}
		printf("\n");
	}	
}
void begin(char *data[5000], int edge_num, char *demand) {
		printf("edge_num:%d\n", edge_num);
		int topo[4];
		int n = 0;
		char *p = demand;
		n = 0;
		int cmdNode;
		while (*p != 0) {
			if (*p == '\r' || *p == '\n') break;
			sscanf(p, "%d", &cmdNode);
			if (n == 0) {
				startNode = cmdNode;
			}
			else
				if (n == 1) {
					endNode = cmdNode;
				}
				else {
					//requireNodes.insert(cmdNode);
					requireNodes[numOfReqNodes++] = cmdNode;					
					isNodesRequire[cmdNode] = true;
				}
				n++;
				p++;
				while (*p >= '0' && *p <= '9') p++;
				if (*p != 0) p++;

		}
	/* init the graph */
	
    for(int i = 0; i < edge_num; i++) {
        p = data[i];		
		n = 0;
		while (*p != 0) {
			if (*p == '\r' || *p == '\n') break;
			sscanf(p, "%d", &topo[n++]);
			p++;
			while (*p >= '0' && *p <= '9') p++;
			if (*p != 0) p++;

		}
		/* 根据wiki: https://en.wikipedia.org/wiki/Travelling_salesman_problem 
		   介绍的方法，将非对称tsp转换为对称tsp求解 */ 
		int linkId = topo[0], from = topo[1], to = topo[2], cost = topo[3];
		/* from相当于原始的A, to相当于原始的B  
		   那么Distance[A'][B] = Distance[B][A'] = Distance[A][B]
		   在这里，取A'=A+600;
		   由于将矩阵扩大了一倍，所以原来的位置依然保存的是原本的数据。
		*/
		int A = from, B = to;
		int A_p = A + 600, B_p = B + 600;
		// 起点不应该保存入度，终点不应该保存出度
		if (cost < matrix[from][to].cost) {
			// TODO:如果cost>原本的，应该不保存！继续continue!!
			matrix[from][to].linkid = linkId;
			matrix[from][to].cost = cost;			
		}else{
			continue;
		}
		/* 将A_p A B_p B的出度保存到node_map_for_ls里
		   A_p是B的出度，B也是A_p的出度
		*/
		//printf("from=%d to=%d\n", from, to);
		if(from == endNode || to == startNode) {
			// 是从终点出发的或者到达起点的边，都不保存！
			continue;
		}
			add_to_node_map_for_ls(from, to, cost);
		// 这里用的是正常的节点信息
			map<unsigned short, NodeInfo>::iterator fromIt = nodeInfoMap.find(from);
			if (fromIt != nodeInfoMap.end())
			{
				// 出度+1
				fromIt->second.outdegree += 1;
				fromIt->second.outdegree_list.push_back(to);
				if (isNodesRequire[to]) {
					fromIt->second.out_include_reqnode_cnt += 1;
				}
				//printf("from:%d found, outdegree++,=%d\n", fromIt->second.nodeId, fromIt->second.outdegree);
			}
			else {
				//printf("from not found, outdegree = 1\n");
				NodeInfo nodeInfo;
				nodeInfo.outdegree = 1;
				nodeInfo.indegree = 0;
				if (isNodesRequire[to]) {
					nodeInfo.out_include_reqnode_cnt = 1;
				}
				else {
					nodeInfo.out_include_reqnode_cnt = 0;
				}
				nodeInfo.outdegree_list.push_back(to);
				nodeInfoMap[from] = nodeInfo;
			}
		
			map<unsigned short, NodeInfo>::iterator toIt = nodeInfoMap.find(to);
			if (toIt != nodeInfoMap.end())
			{
				// 入度+1
				toIt->second.indegree += 1;
				toIt->second.indegree_list.push_back(from);
				//printf("to:%d found, indegree++,=%d\n", toIt->second.nodeId, toIt->second.indegree);
			}
			else {
				//printf("to not found, indegree = 1\n");
				NodeInfo nodeInfo;
				nodeInfo.outdegree = 0;
				nodeInfo.indegree = 1;
				nodeInfo.out_include_reqnode_cnt = 0;
				nodeInfo.indegree_list.push_back(from);
				nodeInfoMap[to] = nodeInfo;
			}
			nodes.insert(from);
			nodes.insert(to);		
		
	}	
	reqNodeMaxOutDegree = nodeMaxIndegree = 0;
	// 统计出入度信息
	for (map<unsigned short, NodeInfo>::iterator it = nodeInfoMap.begin();
	it != nodeInfoMap.end(); it++) {
		unsigned short nodeId = it->first;
		NodeInfo & nodeInfo = it->second;
		if (nodeInfo.indegree == 0) {
			if(nodeId == endNode || isNodesRequire[nodeId]){
				// 无解
				printf("end node:%d indegree is 0\n", nodeId);
				exit(0);
			}
		}
		if(nodeInfo.indegree > nodeMaxIndegree) {
			nodeMaxIndegree = nodeInfo.indegree;
		}
		if (isNodesRequire[nodeId]) {			
			//printf("node=%d indegree=%d outdegree=%d out_include_reqnode=%d\n", nodeId, nodeInfo.indegree,
			//	nodeInfo.outdegree, nodeInfo.out_include_reqnode_cnt);			
			if (nodeInfo.outdegree > reqNodeMaxOutDegree) {
				reqNodeMaxOutDegree = nodeInfo.outdegree;
			}			
			if (nodeInfo.indegree == 1) {
				unsigned short v = nodeInfo.indegree_list[0];
				//	printf("%d has only one indegree:%d\n",nodeId, v);
				mustGoToReqNode[v] = nodeId;
			}
		}
	}
    printf("startNode=%d endNode=%d\n", startNode, endNode);
	printf("nodeMaxIndegree=%d reqNodeMaxOutDegree=%d\n", nodeMaxIndegree, reqNodeMaxOutDegree);
	//printf("nodes size=%d requireNodes size=%d\n", nodes.size(), numOfReqNodes);		
	srand(time(NULL));
	trans_graph_to_atsp();
	initProgram();	
	/* 转换之后就可以用LKH求解了 */
	// key = 在矩阵中从1到n的顺序,value=保存的node_id
	vector<int> node_to_write;
	map<int, int> convert_node_id;
	printf("node to write:\n");
	for(map<unsigned short, NodeInfo>::iterator it = node_map_for_ls.begin(); 
		it != node_map_for_ls.end(); it++)
	{
		unsigned short node = (*it).first;
		printf("%d ", node);
		node_to_write.push_back(node);
	}	
	printf("\n");
	int dimension = node_to_write.size();
	/* 输出图数据  */
	FILE* para_file = fopen("para.tsp", "w");
	if(!para_file){
		printf("open para.atsp error\n");
		//exit(1);
	}
	char *type = "TYPE: TSP\n";
	char *EDGE_WEIGHT_TYPE = "EDGE_WEIGHT_TYPE: EXPLICIT\n";
	char *EDGE_WEIGHT_FORMAT = "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n";
	char *EDGE_WEIGHT_SECTION = "EDGE_WEIGHT_SECTION\n";
	char dimension_str[100];
	sprintf(dimension_str, "DIMENSION: %d\n", dimension);
	fwrite (type , sizeof(char), strlen(type), para_file);
	fwrite (EDGE_WEIGHT_TYPE , sizeof(char), strlen(EDGE_WEIGHT_TYPE), para_file);
	fwrite (dimension_str , sizeof(char), strlen(dimension_str), para_file);
	fwrite (EDGE_WEIGHT_FORMAT , sizeof(char), strlen(EDGE_WEIGHT_FORMAT), para_file);
	fwrite (EDGE_WEIGHT_SECTION , sizeof(char), strlen(EDGE_WEIGHT_SECTION), para_file);
	int nodeid_in_matrix = 1;
	char buffer[50];
	char line_feed[3] = "\n";
	for(int i = 0; i < dimension; i++) {
		// 横向是node_to_write[0]....node_to_write[n]
		// 纵向也是node_to_write[0]....node_to_write[n]
		int from_node = node_to_write[i];
		// 表示在矩阵中第nodeid_in_matrix个点原本是from_node
		convert_node_id[nodeid_in_matrix++] = from_node;
		/* 按node_to_write[0]....node_to_write[n]顺序输出从
		 from_node到node_to_write[j]的权重
		*/
		for(int j = 0; j < dimension; j++) {
			int to_node = node_to_write[j];
			int cost = weightMatrix[from_node][to_node].cost;
			sprintf(buffer, "%d\t", cost);
			fwrite(buffer, sizeof(char), strlen(buffer), para_file);
		}
		fwrite(line_feed, sizeof(char), strlen(line_feed), para_file);
	}
	fclose(para_file);
	
	bool found = construct_first_tour();
	vector<int> result;
	local_search(best_so_far_ant->tour, result, best_so_far_ant->isForward);
	FILE* tour_file = fopen("start_tour.txt", "w");
	if(!tour_file){
		printf("open start_tour.txt error\n");
		//exit(1);
	}
	printf("result:\n");
	for(int i = 0; i < result.size()-1; i++){ 
		printf("%d ", result[i]);
		sprintf(buffer, "%d\t", result[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), tour_file);
	}
	fclose(tour_file);
	printf("\n");
	return;	
			
	found = construct_first_tour();
	//bool found = findPathUsingSCC();	
	printf("found:----------------\n\n");
	printf("length=%d requireCnt=%d path:\n", best_so_far_ant->tour_length, best_so_far_ant->requireCnt);
	for(int j = 0; j < best_so_far_ant->tour.size(); j++) {
			printf("%d ", best_so_far_ant->tour[j]);
	}
	printf("\n\n");
	ignoreWeightMode = false;	
	if(found) {					
		if(max_tries > 5) {
			//max_tries = (max_tries >> 1) + 1;
		}
		//max_tries = 30;
		max_tries = numOfReqNodes;
		if(max_tries <= 30){
		   max_tries = 30;	
		}
		max_tours = 100;
		printf("found tour length:%d\n", best_so_far_ant->tour_length);		
		double p_x = exp(log(0.05)/nodes.size()); 
		trail_min = 1. * (1. - p_x) / (p_x * (double)((8) / 2));
		trail_max = 1. / ( best_so_far_ant->tour_length );
		trail_0 = trail_max;
		trail_min = trail_max * trail_min; 
	}
	int loopTime = 0;
	while ((loopTime == 0) || (
		best_so_far_ant->requireCnt < numOfReqNodes && elapsed_time(VIRTUAL) < max_time))
	{	
			//if(loopTime >= 1) break;
			for (int n_try = 0; elapsed_time(VIRTUAL) < max_time &&  n_try <= max_tries ; n_try++ ) {
				if(loopTime > 0) {
					// 增加迭代次数
					max_tours += 10;
					if(max_tours >= 350) max_tours = 350;
				}
				init_try(n_try%10);

				while ( !termination_condition() ) {	
						//rewardForShortPath();
						construct_solutions();						
						update_statistics();
						pheromone_trail_update();  
						punishAndReward();
						//rewardForDijkstraShortPath();
						search_control_and_statistics();
						iteration++;
						// 奖励计算太多，放在一起重新计算total
						compute_total_information();
				}
				//printf("exit_try:%d\n", n_try);
			}	
			loopTime++;
		}
		if(elapsed_time(VIRTUAL) < max_time) {
			search_with_single_ant(100);		
		}
			printf("\n--------- best=%d\n", best_so_far_ant->tour_length);	
			  int tourSize = best_so_far_ant->tour.size();
			  if(best_so_far_ant->isForward){
                   for (int i = 0; i < tourSize - 1; i++) {
                                                //printf("i=%d\n", i);
                        unsigned short from = best_so_far_ant->tour[i],
                        to = best_so_far_ant->tour[i+1];
                        //printf("%d %d \n", from, to);
                         record_result(matrix[from][to].linkid);
                   }   
               }else{
    
                  for (int i = tourSize-1; i > 0; i--) {
                       //printf("i=%d\n", i);
                        unsigned short from = best_so_far_ant->tour[i],
                         to = best_so_far_ant->tour[i-1];
                        // printf("%d %d \n", from, to);
                         record_result(matrix[from][to].linkid);
                  }   
            }   
	//delete[] isNodesRequire;
	delete[] ant;
	delete best_so_far_ant;
	delete iteration_best_ant;
	delete best_reverse_ant;
	delete sub_opt_so_far_ant;
	delete[] visited_for_ls;

}
