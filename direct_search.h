#ifndef DIRECT_SEARCH_H
#define DIRECT_SEARCH_H

#include "ants.h"
#include <stack>
#include <vector>
#include <stdio.h>
#include <set>
#include <queue>
#include <map>
#include <algorithm>
#include <utility>
#include <math.h>
#include "timer.h"

struct HeapNode {
	unsigned short dist, nodeid;
	bool operator < (const HeapNode& rhs) const {
		return dist > rhs.dist;
	}
};
// ÓÃdijkstraÇó³öµÄÆðµãµ½ËùÓÐ±Ø¾­µãµÄÂ·¾¶£¬±Ø¾­µãÁ½Á½Ö®¼ä×î¶ÌÂ·¾¶£¬±Ø¾­µãµ½ÖÕµã×î¶ÌÂ·¾¶
// ½á¹û
struct DijkShortPathResult{
	// ½á¹û¼¯±£´æµÄÊÇÄÄ¸ö½Úµã³ö·¢µÄµ¥Ôª×î¶ÌÂ·¾¶µÄ½á¹û
	unsigned short nodeId;
	/* key = to_node_id,Ò²¾ÍÊÇµ±Ç°½ÚµãÒªµ½´ïµÄ½ÚµãµÄ×î¶ÌÂ·¾¶ */
	std::map<unsigned short , std::vector<unsigned short> > result;
};



void shortest_path_between_reqnode(std::vector<std::vector<unsigned short> >& all_result);

// snode=Æðµã
void dijkstra(unsigned short snode, bool  *visit, int *dist, unsigned short *parent);


/* ÐÂÔö¡°ÆøÎ¶¡±×÷ÎªÆô·¢Òò×Ó */
static void dfs_smell(const unsigned short& max_range,const unsigned short pre_node, 
			   const unsigned short cur_node, unsigned short cur_range, bool *visited,
			   const unsigned short& req_node);

void init_smell();

void search_with_single_ant(const unsigned short& max_deep);
/* max_deep:×î´óËÑË÷Éî¶È£¬Ò²¾ÍÊÇ×î¶àËÑ³ö¶àÉÙ¸ö½á¹û¾Í½áÊø */
void find_best_with_smell_dfs(const unsigned short& cur_node, const unsigned short& cur_len, unsigned short& max_len,
						  const unsigned short& pre_node, unsigned short* parent, unsigned short req_cnt,
						  bool* visited, std::vector<unsigned short>& result, unsigned short& deep, const unsigned short& max_deep);


void bfs_find_req_cnt(const unsigned short& cur_node, unsigned short& req_cnt, bool* visit);

#endif