#include "stdafx.h"
#include "direct_search.h"

using namespace std;


bool* visit_for_bfs; // for bfs_find_req_cnt

void shortest_path_between_reqnode(vector<vector<unsigned short> >  &result){
	int node_size = nodes.size()+50;
	bool  *visit = new bool[node_size];
	// Æðµãµ½¸÷µãµÄ¾àÀë
	int *dist = new int[node_size];
	// parent[i]=iµÄÉÏÒ»¸ö½Úµã
	unsigned short *parent = new unsigned short[node_size];
	vector<unsigned short> node_to_find;
	for(int i = 0; i < numOfReqNodes; i++){
		node_to_find.push_back(requireNodes[i]);
	}
	node_to_find.push_back(startNode);
	node_to_find.push_back(endNode);
	for(int i = 0; i < node_to_find.size(); i++) {
		dijkstra(node_to_find[i], visit, dist, parent);				
		for(int j = 0; j < node_to_find.size(); j++){
			if(j == i) continue;
			vector<unsigned short> path;
			unsigned short nextNode = node_to_find[j];
			//printf("curNode:%d node_to_find:%d\n", node_to_find[i], nextNode);
			while(nextNode != 900) {
				path.push_back(nextNode);
				//printf("%d ", nextNode);
				nextNode = parent[nextNode];
			}		
			//printf("\n");
			result.push_back(path);
		}
	}
	delete[] dist;
	delete[] parent;

}

void dijkstra(unsigned short snode, bool  *visit, int *dist, unsigned short *parent){
	int node_size = nodes.size();
	for(int i = 0; i < node_size; i++) {
		visit[i] = false;
		dist[i] = 12000;
		parent[i] = 900;
	}
	dist[snode] = 0;
	priority_queue<HeapNode> queue;
	HeapNode head;
	head.dist = 0;
	head.nodeid = snode;
	queue.push(head);
	while(!queue.empty()) {
		HeapNode x = queue.top(); queue.pop();
		unsigned short nodeid = x.nodeid;
		if(visit[nodeid]) continue;
		visit[nodeid] = true;
		const NodeInfo &nodeInfo = nodeInfoMap[nodeid];
		for(int i = 0; i < nodeInfo.outdegree; i++) {
			unsigned short nextNode = nodeInfo.outdegree_list[i];
			int to_dist = matrix[nodeid][nextNode].cost;
			if(dist[nextNode] > dist[nodeid] + to_dist) {
				dist[nextNode] = dist[nodeid] + to_dist;
				parent[nextNode] = nodeid;
				HeapNode headNode;
				headNode.dist = dist[nextNode];
				headNode.nodeid = nextNode;
				queue.push(headNode);
			}
		}
	}
}

/* ³õÊ¼»¯ÆøÎ¶,ÒÔ±Ø¾­µãÎªÒÀ¾Ý£¬±Ø¾­µã¾àÀë30ÒÔÄÚµÄµãµÄÆøÎ¶¶¼±»ÏàÓ¦µÄÔö¼Ó£»
  Ê¹ÓÃdfs´Ó±Ø¾­µã¿ªÊ¼ÍùÈë¶ÈËÑË÷¡£	*/
void init_smell() {
	// ³õÊ¼ÆøÎ¶
	double init_smell = 0.01;
	unsigned short max_range = 25;
	bool *visited = new bool[MAX_NODES];
	//for(int i = 0; i < MAX_EDGES; i++) 
	//	smell[i] = init_smell;		
	
	for(int n = 0; n < numOfReqNodes; n++){
		for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++){
				visited[*it] = false;
		//	for(set<unsigned short>::iterator to = nodes.begin(); to != nodes.end(); to++){
		//		if(matrix[*it][*to].cost == 0) continue;
		//		visited[matrix[*it][*to].linkid] = false;
		//	}
		}			
		unsigned short cur_node = requireNodes[n];
		//printf("reqNode:%d\n", cur_node);
		dfs_smell(max_range, cur_node, cur_node, 0, visited, cur_node);	
		/*	
		for(set<unsigned short>::iterator it = nodes.begin(); it != nodes.end(); it++){
			const map<unsigned short, double>& smell_info = smell[*it];
			const map<unsigned short, double>::const_iterator m_it = smell_info.find(cur_node);
			if(m_it != smell_info.end()) {
				printf("node=%d smell=%f\n", *it, m_it->second);
			}
			
			for(set<unsigned short>::iterator to = nodes.begin(); to != nodes.end(); to++){
				if(*to == *it) continue;
				if(matrix[*it][*to].cost == 0) continue;
				const map<unsigned short, double>& smell_info = smell[matrix[*it][*to].linkid];
				const map<unsigned short, double>::const_iterator m_it = smell_info.find(cur_node);
				if(m_it != smell_info.end()) {
					printf("%d -> %d smell=%f\n", *it, *to, m_it->second);
				}
			}
			
		}
		*/
	}
	delete[] visited;

}

static void dfs_smell(const unsigned short& max_range,const unsigned short pre_node, 
			   const unsigned short cur_node, unsigned short cur_range, bool *visited,
			   const unsigned short& req_node) {
				  // unsigned short link_id = 4805;
				  // if(pre_node != cur_node) {
					//   link_id = matrix[cur_node][pre_node].linkid;
				 //  }
				  // if(link_id != 4805 && visited[link_id]) return;
				   if(visited[cur_node]) return;
				   if(cur_node == startNode || cur_node == endNode || cur_range >= max_range) {
					   // µ½ÁËÆðµã»òÕßÖÕµã¾ÍÍ£Ö¹,»òÕßµ½ÁË×î´ó·¶Î§
					   visited[cur_node] = true;
					   // ÓÉÓÚÊÇÓÃÈë¶ÈÀ´ËãµÄ£¬ËùÒÔ·½ÏòÆäÊµÊÇcur_node->pre_node
					   int cost = matrix[cur_node][pre_node].cost;
					  double weight = 1.0 / (cur_range+0.0001);
					   map<unsigned short, double>& smell_info = smell[cur_node];
					   if(smell_info.find(req_node) != smell_info.end()) {
						   // ÒÑ¾­´æÔÚÁË
						   smell_info[req_node] += weight;
					   }else{
						   smell_info[req_node] = weight;
					   }					  		   
					   return;
				   }				  
				   visited[cur_node] = true;
				   const NodeInfo& node_info = nodeInfoMap[cur_node];
				   vector<unsigned short>::const_iterator it = node_info.indegree_list.begin();
				   for(; it != node_info.indegree_list.end(); it++) {
					   unsigned short next_node = *it;
					   //unsigned short local_link_id = matrix[next_node][cur_node].linkid;
					   if(visited[next_node] == false) {
						   int cost = matrix[next_node][cur_node].cost;
						   if(cur_range+cost <= 20) 
							dfs_smell(max_range, cur_node, next_node, cur_range+cost, visited, req_node);
					   }
				   }
				   if(cur_node != pre_node) {
					   // ÓÉÓÚÊÇÓÃÈë¶ÈÀ´ËãµÄ£¬ËùÒÔ·½ÏòÆäÊµÊÇcur_node->pre_node
					   int cost = matrix[cur_node][pre_node].cost;
					  // double weight = 1.0 / (cost + 0.0001);
					   double weight = 1.0 / (cur_range + cost + 0.0001);
					   map<unsigned short, double>& smell_info = smell[cur_node];
					   if(smell_info.find(req_node) != smell_info.end()) {
						   // ÒÑ¾­´æÔÚÁË
						   smell_info[req_node] += weight;
					   }else{
						   smell_info[req_node] = weight;
					   }		
				   }

}

// Ö»ÓÃÒ»Ö»ÂìÒÏËÑË÷(ÆäÊµ¾ÍÊÇÓÃdfs+»ØËÝ+ÆøÎ¶À´³¢ÊÔËÑ×îÓÅ½â£©
void search_with_single_ant(const unsigned short& max_deep) {
	// È¨ÖØµÄÉÏ½ç
	unsigned short max_len = best_so_far_ant->tour_length;
	printf("best_so_far_ant->tour_length=%d\n", best_so_far_ant->tour_length);
	init_pheromone_trails(trail_min);
	global_acs_pheromone_update(best_so_far_ant);
	compute_total_information();
	int max_size = nodes.size() + 10;
	bool* visited = new bool[max_size];
	
	unsigned short *parent = new unsigned short[max_size];
	vector<unsigned short> result;
	for(int j = 0; j < max_size; j++) visited[j] = false;			
	for(int j = 0; j < max_size; j++) parent[j] = 900;
	// TODO£º¿ÉÒÔ³¢ÊÔÓÃÕÒ³öÀ´µÄ½â¶ÔÐÅÏ¢ËØ½øÐÐÇ¿»¯£¬È»ºó
	//  ¸ù¾ÝÐÅÏ¢ËØ+ÆøÎ¶¹²Í¬¾ö¶¨ËÑÄ³¸öµãµÄ¸ÅÂÊ 
	unsigned short deep = 0;
	find_best_with_smell_dfs(startNode, 0, max_len, 900, parent, 0, visited, result, deep, max_deep);
	printf("max_len=%d path size=%d\n", max_len, result.size());
	if(max_len < best_so_far_ant->tour_length) {
		best_so_far_ant->tour.clear();
		best_so_far_ant->isForward = true;
		best_so_far_ant->curNode = startNode;
		best_so_far_ant->tour_length = max_len;
		for(int i = result.size()-1; i >= 0 ; i--) {
		//	printf("%d ", result[i]);
			best_so_far_ant->tour.push_back(result[i]);
		}
	}	
	//printf("\n");
	delete[] visited;
	delete[] parent;
}
// µÚÒ»Î¬±£´æµÄÊÇnodeid£¬µÚ¶þÎ¬ÊÇ¶ÔÓ¦µÄsmell
typedef pair<unsigned short, double> pud;
bool cmp(const pud &a, const pud&b)
{
	if (a.second > b.second)
        return true;
    else return false;
}


void find_best_with_smell_bfs()
{
	
}

/* ÓÃdfsËÑË÷×î¼ÑÂ·¾¶ */
void find_best_with_smell_dfs(const unsigned short& cur_node, const unsigned short& cur_len, unsigned short& max_len,
						  const unsigned short& pre_node, unsigned short* parent, unsigned short req_cnt,
						  bool* visited, vector<unsigned short>& result, unsigned short& deep, const unsigned short& max_deep) 
{
	if(elapsed_time(VIRTUAL) >= 9) return;	
	//if(deep > max_deep) return;
	if(visited[cur_node]) return;
	if(req_cnt >= numOfReqNodes && cur_node == endNode) {
		//printf("cur_len=%d\n", cur_len);
		deep++;
		// µ½´ïÖÕµãÁË
		if(cur_len < max_len) {
			max_len = cur_len;
			result.clear();
			unsigned short to_p = parent[cur_node];
			result.push_back(cur_node);
			// Èç¹ûto_nodeµÄparent != 900,ËµÃ÷µ½´ïÁË£¬ËµÃ÷ÕÒµ½ÁËÒ»Ìõ¸ü¶ÌµÄÂ·¾¶		
			while(to_p != 900) {
				result.push_back(to_p);
				to_p = parent[to_p];
			}
		}
		return;
	}
	visited[cur_node] = true;
	const NodeInfo& node = nodeInfoMap[cur_node];
	int nn_ants = node.outdegree;
	pud  prob_ptr[9];	
	vector<unsigned short>::const_iterator adjlist_begin = node.outdegree_list.begin();
	vector<unsigned short>::const_iterator adjlist_end = node.outdegree_list.end();
	int idx = 0;
	double   rnd, partial_sum = 0., sum_prob = 0.0;
	// ¼ÆËã¸ÅÂÊ
	for(; adjlist_begin != adjlist_end; adjlist_begin++) 
	{
		unsigned short to = *adjlist_begin;
		int cost = matrix[cur_node][to].cost;						
		if(visited[to] || cur_node == to || (cost+cur_len) >= max_len) {
			// ÒÑ¾­·ÃÎÊ¹ý¡¢ÊÇµ±Ç°×Ô¼º¡¢È¨ÖØÒÑ¾­³¬³öÉÏÏÞ
			//prob_ptr[idx++] = make_pair(to, 0);
		}				
		else{
			/*
			bool is_node_ok = true;
			if(req_cnt > 1) {
				// ¶ÔÃ¿Ò»¸öÏàÁÚµÄµã£¬ÓÃbfsËÑË÷¿´´ÓÕâ¸öµã³ö·¢£¬
				//	ÊÇ·ñÄÜ¾­¹ýËùÓÐ±Ø¾­µã£¬Èç¹û²»ÐÐ£¬ÄÇ¾Í²»×ßÕâ¸öµã
				unsigned short tmp_req_cnt = req_cnt;
				int max_size = nodes.size() + 5;
				for(int  n = 0; n < max_size; n++) {
					visit_for_bfs[n] = visited[n];
				}
				bfs_find_req_cnt(to, tmp_req_cnt, visit_for_bfs);
				if(tmp_req_cnt + req_cnt < numOfReqNodes)
					is_node_ok = false;
			}
			if(is_node_ok == false) continue;
			*/
			// TODO£º¶ÔÃ¿Ö»ÂìÒÏÀ´Ëµ£¬¼ÓÈëµÄÆøÎ¶ÊÇ²»Ò»ÑùµÄ
			double smell_weight = 0.0001;
			bool has_smell = false;
			bool in_best_tour = best_so_far_ant->visited[to];
			const map<unsigned short, double>& smell_info = smell[to];
			map<unsigned short, double>::const_iterator m_it = smell_info.begin();
			for(; m_it != smell_info.end(); m_it++) {
				unsigned short smell_req_node = m_it->first;
				if(isNodesRequire[smell_req_node] && visited[smell_req_node] == false) {
					//smell_weight += m_it->second;
					has_smell = true;
					// find the larger one
					if(m_it->second > smell_weight)
						smell_weight = m_it->second;
				}
			}
			//prob_ptr[idx] = total[cur_node][to];	
			
			double weight = pow(smell_weight, gama)*total[cur_node][to];
			//double weight = pow(smell_weight, gama);
			prob_ptr[idx] = make_pair(to, weight);				
			idx++;
		}
	}
	sort(prob_ptr, prob_ptr + nn_ants, cmp);
	if(idx >= 4) idx = 4;
	for(int i = 0; i < idx; i++) {
		// Ñ¡ÔñÏÂÒ»¸öµã			
		unsigned short nextNode = prob_ptr[i].first;
		unsigned short is_req = 0;
		if(isNodesRequire[nextNode]) {
			is_req++;
		}		
		int total_len = cur_len + matrix[cur_node][nextNode].cost;
		parent[nextNode] = cur_node;
		find_best_with_smell_dfs(nextNode, total_len,max_len, cur_node, parent, 
							req_cnt+is_req, visited, result, deep, max_deep);
		visited[nextNode] = false;
	}
}

/* Í¨¹ýbfsÀ´±éÀúÍ¼£¬À´Í³¼Æ´Ós_node³ö·¢£¬ÔÚvisitÊý×éµÄÇé¿öÏÂ£¬
	ÄÜÕÒµ½¶àÉÙ¸ö±Ø¾­µã¡£*/
void bfs_find_req_cnt(const unsigned short& cur_node, unsigned short& req_cnt, bool* visit)
{	
	queue<unsigned short> node_q;
	node_q.push(cur_node);
	while(node_q.empty() == false) {
		if(req_cnt >= numOfReqNodes) break;
		unsigned short next_node = node_q.front();
		node_q.pop();
		if(visit[next_node]) continue;		
		visit[next_node] = true;
		if(isNodesRequire[next_node]) req_cnt++;
		const NodeInfo& node = nodeInfoMap[next_node];
		for(int i = 0; i < node.outdegree_list.size(); i++) {
			unsigned short to_node = node.outdegree_list[i];
			if(visit[to_node] == false){
				node_q.push(to_node);
			}
		}
	}
}
