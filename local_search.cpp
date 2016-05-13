#include "stdafx.h"
#include "local_search.h"
#include <assert.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <algorithm>

using namespace std;


int ls_flag;          /* indicates whether and which local search is used */ 
int nn_ls;            /* maximal depth of nearest neighbour lists used in the 
			      local search */ 
bool dlb_flag = true;  /* flag indicating whether don't look bits are used. I recommend 
			      to always use it if local search is applied */

unsigned short *pos;               /* positions of cities in tour */ 
bool *dlb;               /* vector containing don't look bits */ 
unsigned short *h_tour;            /* help vector for performing exchange move */ 
unsigned short *hh_tour;           /* help vector for performing exchange move */ 
 unsigned short *random_vector;

/* 初始化local search会用到的数组 */
void init_ls() {
	int size = MAX_NODES*4;
	pos = new unsigned short[size];
    dlb = new bool[size];
    h_tour = new unsigned short[size];
    hh_tour =new unsigned short[size];
	random_vector = new unsigned short[size];  
}

void free_ls() {
	delete[] h_tour;
	delete[] hh_tour;
	delete[] pos;
	delete[] dlb;  
	delete[] random_vector;
}
void generate_random_permutation( int n )
{
	vector<unsigned short> number;
	for(int i = 0; i < n; i++) {
		number.push_back(i);
	}
	random_shuffle(number.begin(), number.end());	
	for(int i = 0; i <number.size(); i++) random_vector[i] = number[i];
}

/**
**	visited:用来表明某个点是否是在当前路径上的
	n表示tour的点的个数
  ghost_cnt表示不在路径上的点，=图的所有结点nodes.size()-在路径上的点
*/
//void three_opt_first( vector<unsigned short>& tour, const unsigned short& n, const bool* visited,
//					 const unsigned short& ghost_cnt)
void three_opt_first( vector<unsigned short>& tour, const unsigned short& n)

/*    
      FUNCTION:       3-opt the tour
      INPUT:          pointer to the tour that is to optimize
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 3-opt
      COMMENT:        this is certainly not the best possible implementation of a 3-opt 
                      local search algorithm. In addition, it is very lengthy; the main 
		      reason herefore is that awkward way of making an exchange, where 
		      it is tried to copy only the shortest possible part of a tour. 
		      Whoever improves the code regarding speed or solution quality, please 
		      drop me the code at stuetzle no@spam ulb.ac.be
		      The neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers available online from David S. Johnson.
*/
{
    /* In case a 2-opt move should be performed, we only need to store opt2_move = true,
       as h1, .. h4 are used in such a way that they store the indices of the correct move */

    unsigned short   c1, c2, c3;           /* cities considered for an exchange */
    unsigned short   s_c1, s_c2, s_c3;     /* successors of these cities        */
    unsigned short   p_c1, p_c2, p_c3;     /* predecessors of these cities      */   
    unsigned short   pos_c1, pos_c2, pos_c3;     /* positions of cities c1, c2, c3    */
    unsigned short   i, j, h, g, l;
    int   help;
	bool  improvement_flag;
    unsigned short   h1=0, h2=0, h3=0, h4=0, h5=0, h6=0; /* memorize cities involved in a move */
    int  diffs, diffp;
    bool   between = false; 
    unsigned short   opt2_flag;  /* = true: perform 2-opt move, otherwise none or 3-opt move */
    unsigned short   move_flag;  /* 
			      move_flag = 0 --> no 3-opt move 
			      move_flag = 1 --> between_move (c3 between c1 and c2)
			      move_flag = 2 --> not_between with successors of c2 and c3
			      move_flag = 3 --> not_between with predecessors of c2 and c3
			      move_flag = 4 --> cyclic move 
			   */
    int gain, move_value, radius, add1, add2;
    int decrease_breaks;    /* Stores decrease by breaking two edges (a,b) (c,d) */
    int val[3];
    unsigned short n1, n2, n3;

	//unsigned short pos[610], h_tour[610], hh_tour[610];
	//bool dlb[610];
   

    for ( i = 0 ; i < n ; i++ ) {
	    // 保存tour[i]的点的下标
		pos[tour[i]] = i;
		dlb[tour[i]] = false;
    }
	//printf("pos:\n");
	//for ( i = 0 ; i < n ; i++ ) {
	//   printf("node:%d pos:%d\n", tour[i], pos[tour[i]]);
   // }
//	dlb[endNode] = true;
//	dlb[endNode + 600] = true;
    improvement_flag = true;
	generate_random_permutation(n);
	//printf("random_vector:\n");
	//for(int i = 0; i < n; i++) 
	//	printf("%d ", random_vector[i]);
	//printf("\n\n");
    while ( improvement_flag ) {
	move_value = 0;
	improvement_flag = false;
	//dlb_flag = false;
	for ( l = 0 ; l < n ; l++ ) {
		// 随机选一个点作为c1
		unsigned short random_idx = random_vector[l];
		assert(random_idx < tour.size() && random_idx >= 0);
	    c1 = tour[random_idx];		
		// >= VIRT_OFFSET的都是虚拟点。虚拟点不要作为c1
		if ((dlb_flag && dlb[c1]))
			continue;
	    opt2_flag = false;

	    move_flag = 0;		// no 3-opt move 
	    pos_c1 = pos[c1];
		if(pos_c1 >= tour.size()) {
				printf("c1=%d c2=%d pos_c1=%d\n", c1, c2, pos_c1);
		}
	    // 路径tour里，c1的下一个节点
		int pos_c1_p1 = pos_c1+1;
		if(pos_c1_p1 >= tour.size()) {
			pos_c1_p1 = 0;
		}
	    s_c1 = tour[pos_c1_p1];
	    if (pos_c1 > 0)
			 // 路径tour里，c1的上一个节点
			p_c1 = tour[pos_c1-1];
	    else 
			p_c1 = tour[n-1];

	    h = 0;    /* Search for one of the h-nearest neighbours */
		nn_ls = node_map_for_ls[c1].outdegree;
	    while ( h < nn_ls ) {
			// 这个应该是说是c1的邻接表里面的点，而不像p_c1是tour里的下一个点
			// 注意 使用的是为local_search定制的node_map_for_ls
			assert(h < node_map_for_ls[c1].outdegree_list.size());
			c2   = node_map_for_ls[c1].outdegree_list[h];  /* second city, determine its position */						
			pos_c2 = pos[c2];
			 // 路径tour里，c1的下一个节点
			int pos_c2_p1 = pos_c2+1;
			if(pos_c2_p1 >= tour.size()) {
				pos_c2_p1 = 0;
			}
			// 但由于普通的tsp问题是把全部点都遍历了，所以不管是哪个点，都会在tour里面存在			
			s_c2 = tour[pos_c2_p1];
			if (pos_c2 > 0)
				p_c2 = tour[pos_c2-1];
			else 
				p_c2 = tour[n-1];
			
			diffs = 0; diffp = 0;
			// 以radius为半径搜索c1的最近邻点（距离超过radius的邻接点不会考虑）
			radius = weightMatrix[c1][s_c1].cost;
			add1   = weightMatrix[c1][c2].cost;
			int dist_c2_sc2 = weightMatrix[c2][s_c2].cost;
			int dist_c2_pc2 = weightMatrix[c2][p_c2].cost;			
			/* Here a fixed radius neighbour search is performed */
			if (radius != 10000 && radius > add1 ) {
			//if (radius > add1 ) {
				// 将(c1,s_c1)(c2,s_c2)这两条边删掉之后，decrease了多少
				//decrease_breaks = - radius - weightMatrix[c2][s_c2].cost;
				decrease_breaks = - radius - dist_c2_sc2;
				/*(c1,s_c1)(c2,s_c2)break掉之后，将(c1,c2)(s_c1.s_c2)连接起来，看diffs是多少
				 如果diffs是负数，说明这种重连方式有效*/
				diffs =  decrease_breaks + add1 + weightMatrix[s_c1][s_c2].cost;
				/* 这里是另一种break的方式:将(c1,s_c1)(c2,p_c2)删掉，然后连接(c1,p_c2)
				( s_c1,c2),然后再求出这种连接方式的增益 */
				diffp =  - radius - dist_c2_pc2 + 
				weightMatrix[c1][p_c2].cost + weightMatrix[s_c1][c2].cost;
			}
			else {
				// TODO;如果邻接表是按升序排列的，则可以break，如果无序应该continue;
				//h++;
				//continue;		
				break;
			}
			if ( p_c2 == c1 )  /* in case p_c2 == c1 no exchange is possible */
				diffp = 0;
			
			// 这些都是2-opt!!不要赋值！！忽略掉！
			if ( (diffs < move_value) || (diffp < move_value) ) {
				improvement_flag = true; 
				if (diffs <= diffp) { 
					// (c1,c2)(s_c1.s_c2)连接的效果更好
					// h1和h3相连，h2和h4相连
					// 这里的意思只是说如果是2-opt交换，那么就用
					// 这些边交换！！
					h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
					move_value = diffs; 
					opt2_flag = true; 
					//opt2_flag = false;
					move_flag = 0;
					//    	    goto exchange; 
				} else {
					// h1和h3相连，h2和h4相连
					// 只对2-opt适用而已！
					h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; 
					move_value = diffp;  
					opt2_flag = true; 
					//opt2_flag = false;
					move_flag = 0;
					//    	    goto exchange; 
				}
			}
			
			/* Now perform the innermost search */
			g = 0;
			// 找s_c1的邻接点
			int g_nn_ls = node_map_for_ls[s_c1].outdegree;
			while (g < g_nn_ls) {
				assert(g < node_map_for_ls[s_c1].outdegree_list.size());
				c3   = node_map_for_ls[s_c1].outdegree_list[g];						
				pos_c3 = pos[c3];
				if(pos_c3+1 >= tour.size()) {
					printf("c1=%d c2=%d c3=%d pos_c3=%d\n", c1, c2, c3, pos_c3);
				}
				int pos_c3_p1 = pos_c3+1;
				if(pos_c3_p1 >= tour.size()) {
					pos_c3_p1 = 0;
				}				
				s_c3 = tour[pos_c3_p1];
				if (pos_c3 > 0)
					p_c3 = tour[pos_c3-1];
				else 
					p_c3 = tour[n-1];
		  
				if ( c3 == c1 ) {
					g++;
					continue;
				}
				else {
					add2 = weightMatrix[s_c1][c3].cost;
					/* Perform fixed radius neighbour search for innermost search */
					/*  decrease_breaks + add1 < add2  这个判断条件超级不合理！！
						这个条件相当于只允许break掉(c2,s_c2)来判断收益，但是在这个if条件里面，
						却又允许break掉的是(c2,p_c2)！！如果break掉的是(c1,s_c1), (c2,p_c2), (c3,s_c3) and 
				   add edges (c1,c2), (c3,s_c1), (p_c2,s_c3) 是有效果的话，这个条件就会导致考虑不到
				   这种情况！！ */
					//if ( decrease_breaks + add1 < add2 ) {
					if (1) {
			  
						if ( pos_c2 > pos_c1 ) {
							if ( pos_c3 <= pos_c2 && pos_c3 > pos_c1 )
								/* c3是在c2和c1之间 */
								between = true;
							else 
								between = false;
						}
						else if ( pos_c2 < pos_c1 )
							if ( pos_c3 > pos_c1 || pos_c3 < pos_c2 )
								between = true;
							else 
								between = false;
						else {
							printf(" Strange !!, pos_1 %ld == pos_2 %ld, \n",pos_c1,pos_c2);
						}
						int dist_c3_sc3 = weightMatrix[c3][s_c3].cost;
						int dist_c3_pc3 = weightMatrix[c3][p_c3].cost;
						int dist_pc3_c3 = weightMatrix[p_c3][c3].cost;						
						if ( between ) {
							/* We have to add edges (c1,c2), (c3,s_c1), (p_c3,s_c2) to get 
							   valid tour; it's the only possibility */
				      
							gain = decrease_breaks - dist_c3_pc3 +
								add1 + add2 + weightMatrix[p_c3][s_c2].cost;
				      
							/* check for improvement by move */
							if ( gain < move_value ) {
								// 这里才是3-opt交换时，要交换的顺序！
								improvement_flag = true; /* g = neigh_ls + 1; */
								move_value = gain;
								/* 如果3-opt交换有效果，那么就不要2-opt交换！
								   如果gain>more_value，说明用3-opt交换没有效果，那么
								   进行的是2-opt交换！
								*/
								opt2_flag = false;
								move_flag = 1;
								/* store nodes involved in move */
								h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = p_c3; h6 = c3;
								goto exchange;
							} 
						}
						else {   /* not between(pos_c1,pos_c2,pos_c3) */
			      
						/* We have to add edges (c1,c2), (s_c1,c3), (s_c2,s_c3) */
						// 将(c1,s_c1), （c2, s_c2)，(c3,s_c3)删掉
						// 用(c1,c2),(s_c1,c3),(s_c2,s_c3)代替
						gain = decrease_breaks - dist_c3_sc3 +
							add1 + add2 + 
							weightMatrix[s_c2][s_c3].cost;
			      
						if ( pos_c2 == pos_c3 ) {
							gain = 20000;
						}
			      
						/* check for improvement by move */
						if ( gain < move_value ) {
							improvement_flag = true; /* g = neigh_ls + 1; */
							move_value = gain;
							opt2_flag = false;
							 // move_flag = 2 --> not_between with successors of c2 and c3
							move_flag = 2;
							/* store nodes involved in move */
							// 连接方式是:(h1,h3),(h2,h5),(h4,h6)
							h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = c3; h6 = s_c3;
							goto exchange;
						}
			      
						/* or add edges (c1,c2), (s_c1,c3), (p_c2,p_c3) */
						gain = - radius - weightMatrix[p_c2][c2].cost 
							- dist_pc3_c3 +
							add1 + add2 +  weightMatrix[p_c2][p_c3].cost;
			      
						if ( c3 == c2 || c2 == c1 || c1 == c3 || p_c2 == c1 ) {
							gain = 2000000;
						}
			      
						if ( gain < move_value ) {
							improvement_flag = true;
							move_value = gain;
							opt2_flag = false;
							move_flag = 3;
							h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = p_c3; h6 = c3;
							goto exchange;
						}
			      
						/* Or perform the 3-opt move where no subtour inversion is necessary 
						   i.e. delete edges (c1,s_c1), (c2,p_c2), (c3,s_c3) and 
						   add edges (c1,c2), (c3,s_c1), (p_c2,s_c3) */
			      
						gain = - radius - weightMatrix[p_c2][c2].cost - 
							dist_c3_sc3 + add1 + add2 + weightMatrix[p_c2][s_c3].cost;
			      
						/* check for improvement */
						if ( gain < move_value ) {
							improvement_flag = true;
							move_value = gain;
							opt2_flag = false;
							move_flag = 4;
							improvement_flag = true;
							/* store nodes involved in move */
							h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = c3; h6 = s_c3; 
							goto exchange;
						}
						}
					}
					else{
						//g++;
					}
				}
				g++;
			}
			h++;
	    }
	    if ( move_flag || opt2_flag ) {
			exchange:
			move_value = 0;

			/* Now make the exchange */
			if ( move_flag ) {
				dlb[h1] = false; dlb[h2] = false; dlb[h3] = false; 
				dlb[h4] = false; dlb[h5] = false; dlb[h6] = false;
				pos_c1 = pos[h1]; pos_c2 = pos[h3]; pos_c3 = pos[h5];
		  
				if ( move_flag == 4 ) {

					if ( pos_c2 > pos_c1 ) 
						n1 = pos_c2 - pos_c1;
					else
						n1 = n - (pos_c1 - pos_c2);
					if ( pos_c3 > pos_c2 ) 
						n2 = pos_c3 - pos_c2;
					else
						n2 = n - (pos_c2 - pos_c3);
					if ( pos_c1 > pos_c3 ) 
						n3 = pos_c1 - pos_c3;
					else
						n3 = n - (pos_c3 - pos_c1);
		      
					/* n1: length h2 - h3, n2: length h4 - h5, n3: length h6 - h1 */
					val[0] = n1; val[1] = n2; val[2] = n3; 
					/* Now order the partial tours */
					h = 0;
					help = INT_MIN;
					for ( g = 0; g <= 2; g++) {
						if ( help < val[g] ) {
						help = val[g];
						h = g;
						}
					}
		      
					/* order partial tours according length */
					if ( h == 0 ) {
						/* copy part from pos[h4] to pos[h5]
						   direkt kopiert: Teil von pos[h6] to pos[h1], it
						   remains the part from pos[h2] to pos[h3] */
						j = pos[h4];
						h = pos[h5];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h) {
						i++;
						j++;
						if ( j  >= n )
							j = 0;
						h_tour[i] = tour[j];
						n1++;
						}
			  
						/* First copy partial tour 3 in new position */
						j = pos[h4];
						i = pos[h6];
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						while ( i != pos_c1) {
						i++;
						if ( i >= n )
							i = 0;
						j++;
						if ( j >= n )
							j = 0;
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						}
			  
						/* Now copy stored part from h_tour */
						j++;
						if ( j >= n )
						j = 0;
						for ( i = 0; i<n1 ; i++ ) {
						tour[j] = h_tour[i];
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						tour[n] = tour[0];
					}
					else if ( h == 1 ) {
			  
						/* copy part from pos[h6] to pos[h1]
						   direkt kopiert: Teil von pos[h2] to pos[h3], it
						   remains the part from pos[h4] to pos[h5] */
						j = pos[h6];
						h = pos[h1];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h) {
						i++;
						j++;
						if ( j  >= n )
							j = 0;
						h_tour[i] = tour[j];
						n1++;
						}
			  
						/* First copy partial tour 3 in new position */
						j = pos[h6];
						i = pos[h2];
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						while ( i != pos_c2) {
						i++;
						if ( i >= n )
							i = 0;
						j++;
						if ( j >= n )
							j = 0;
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						}
			  
						/* Now copy stored part from h_tour */
						j++;
						if ( j >= n )
						j = 0;
						for ( i = 0; i<n1 ; i++ ) {
						tour[j] = h_tour[i];
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						tour[n] = tour[0];
					}
					else if ( h == 2 ) {
						/* copy part from pos[h2] to pos[h3]
						   direkt kopiert: Teil von pos[h4] to pos[h5], it
						   remains the part from pos[h6] to pos[h1] */
						j = pos[h2];
						h = pos[h3];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h) {
						i++;
						j++;
						if ( j  >= n )
							j = 0;
						h_tour[i] = tour[j];
						n1++;
						}
	      
						/* First copy partial tour 3 in new position */
						j = pos[h2];
						i = pos[h4];
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						while ( i != pos_c3) {
						i++;
						if ( i >= n )
							i = 0;
						j++;
						if ( j >= n )
							j = 0;
						tour[j] = tour[i];
						pos[tour[i]] = j; 
						}
			  
						/* Now copy stored part from h_tour */
						j++;
						if ( j >= n )
						j = 0;
						for ( i = 0; i<n1 ; i++ ) {
						tour[j] = h_tour[i]; 
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						tour[n] = tour[0];
					}
				}
				else if ( move_flag == 1 ) {
		      
				if ( pos_c3 < pos_c2 ) 
					n1 = pos_c2 - pos_c3;
				else
					n1 = n - (pos_c3 - pos_c2);
				if ( pos_c3 > pos_c1 ) 
					n2 = pos_c3 - pos_c1 + 1;
				else
					n2 = n - (pos_c1 - pos_c3 + 1);
				if ( pos_c2 > pos_c1 ) 
					n3 = n - (pos_c2 - pos_c1 + 1);
				else
					n3 = pos_c1 - pos_c2 + 1;
		      
				/* n1: length h6 - h3, n2: length h5 - h2, n2: length h1 - h3 */
				val[0] = n1; val[1] = n2; val[2] = n3; 
				/* Now order the partial tours */
				h = 0;
				help = INT_MIN;
				for ( g = 0; g <= 2; g++) {
					if ( help < val[g] ) {
						help = val[g];
						h = g;
					}
				}
				/* order partial tours according length */
		      
				if ( h == 0 ) {
			  
					/* copy part from pos[h5] to pos[h2]
					   (inverted) and from pos[h4] to pos[h1] (inverted)
					   it remains the part from pos[h6] to pos[h3] */
					j = pos[h5];
					h = pos[h2];
					i = 0;
					h_tour[i] = tour[j];
					n1 = 1;
					while ( j != h ) {
						i++;
						j--;
						if ( j < 0 )
							j = n-1;
						h_tour[i] = tour[j];
						n1++;
					}
			  
					j = pos[h1];
					h = pos[h4];
					i = 0;
					hh_tour[i] = tour[j];
					n2 = 1;
					while ( j != h) {
						i++;
						j--;
						if ( j < 0 )
							j = n-1;
						hh_tour[i] = tour[j];
						n2++;
					}
			  
					j = pos[h4];
					for ( i = 0; i< n2 ; i++ ) {
						tour[j] = hh_tour[i];
						pos[hh_tour[i]] = j; 
						j++;
						if (j >= n)
							j = 0;
					}
			  
					/* Now copy stored part from h_tour */
					for ( i = 0; i< n1 ; i++ ) {
						tour[j] = h_tour[i]; 
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
					}
					tour[n] = tour[0];
				}
				else if ( h == 1 ) {
			  
					/* copy part from h3 to h6 (wird inverted) erstellen : */
					j = pos[h3];
					h = pos[h6];
					i = 0;
					h_tour[i] = tour[j];
					n1 = 1;
					while ( j != h) {
						i++;
						j--;
						if ( j  < 0 )
							j = n-1;
						h_tour[i] = tour[j];
						n1++;
					}
			  
					j = pos[h6];
					i = pos[h4];
			  
					tour[j] = tour[i];
					pos[tour[i]] = j; 
					while ( i != pos_c1) {
						i++;
						j++;
						if ( j >= n)
							j = 0;
						if ( i >= n)
							i = 0;
						tour[j] = tour[i];
						pos[tour[i]] = j; 
					}
			  
					/* Now copy stored part from h_tour */
					j++;
					if ( j >= n )
						j = 0;
					i = 0;
					tour[j] = h_tour[i];
					pos[h_tour[i]] = j; 
					while ( j != pos_c1 ) {
						j++;
						if ( j >= n )
							j = 0;
						i++;
						tour[j] = h_tour[i];
						pos[h_tour[i]] = j; 
					}
					tour[n] = tour[0];
				}
		      
				else if ( h == 2 ) {			  
					/* copy part from pos[h2] to pos[h5] and
					   from pos[h3] to pos[h6] (inverted), it
					   remains the part from pos[h4] to pos[h1] */
					j = pos[h2];
					h = pos[h5];
					i = 0;
					h_tour[i] =  tour[j];
					n1 = 1;
					while ( j != h ) {
						i++;
						j++;
						if ( j >= n )
							j = 0;
						h_tour[i] = tour[j];
						n1++;
					}
					j = pos_c2;
					h = pos[h6];
					i = 0;
					hh_tour[i] = tour[j];
					n2 = 1;
					while ( j != h) {
						i++;
						j--;
						if ( j < 0 )
							j = n-1;
						hh_tour[i] = tour[j];
						n2++;
					}
			  
					j = pos[h2];
					for ( i = 0; i< n2 ; i++ ) {
						tour[j] = hh_tour[i];
						pos[hh_tour[i]] = j; 
						j++;
						if ( j >= n)
							j = 0;
					}
			  
					/* Now copy stored part from h_tour */
					for ( i = 0; i< n1 ; i++ ) {
						tour[j] = h_tour[i];
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
					}
					tour[n] = tour[0];
					}
				}
				else if ( move_flag == 2 ) {
		      
					if ( pos_c3 < pos_c1 ) 
						// n1保存的是在tour里面从c3到c1有多少个点？
						n1 = pos_c1 - pos_c3;
					else
						n1 = n - (pos_c3 - pos_c1);
					if ( pos_c3 > pos_c2 ) 
						// 在tour里面从c2到c3有多少个点
						n2 = pos_c3 - pos_c2;
					else
						n2 = n - (pos_c2 - pos_c3);
					if ( pos_c2 > pos_c1 ) 
						// 在tour里面从c1到c2有多少个点
						n3 = pos_c2 - pos_c1;
					else
						n3 = n - (pos_c1 - pos_c2);
		      
					val[0] = n1; val[1] = n2; val[2] = n3; 
					/* Determine which is the longest part */
					h = 0;
					help = INT_MIN;
					for ( g = 0; g <= 2; g++) {
						if ( help < val[g] ) {
							help = val[g];
							h = g;
						}
					}
					/* order partial tours according length */
		      
					if ( h == 0 ) {
			  
						/* copy part from pos[h3] to pos[h2]
						   (inverted) and from pos[h5] to pos[h4], it
						   remains the part from pos[h6] to pos[h1] */
						j = pos[h3];
						h = pos[h2];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h ) {
						i++;
						j--;
						if ( j < 0 )
							j = n-1;
						h_tour[i] = tour[j];
						n1++;
						}
			  
						j = pos[h5];
						h = pos[h4];
						i = 0;
						hh_tour[i] = tour[j];
						n2 = 1;
						while ( j != h ) {
						i++;
						j--;
						if ( j < 0 )
							j = n-1;
						hh_tour[i] = tour[j];
						n2++;
						}
			  
						j = pos[h2];
						for ( i = 0; i<n1 ; i++ ) {
						tour[j] = h_tour[i]; 
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
	      
						for ( i = 0; i < n2 ; i++ ) {
						tour[j] = hh_tour[i];
						pos[hh_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						tour[n] = tour[0];
						/*  	      getchar(); */
					}
					else if ( h == 1 ) {
			  
						/* copy part from pos[h2] to pos[h3] and
						   from pos[h1] to pos[h6] (inverted), it
						   remains the part from pos[h4] to pos[h5] */
						j = pos[h2];
						h = pos[h3];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h ) {
						i++;
						j++;
						if ( j >= n  )
							j = 0;
						h_tour[i] = tour[j];
						n1++;
						}
			  
						j = pos[h1];
						h = pos[h6];
						i = 0;
						hh_tour[i] = tour[j];
						n2 = 1;
						while ( j != h ) {
						i++;
						j--;
						if ( j < 0 )
							j =  n-1;
						hh_tour[i] = tour[j];
						n2++;
						}
						j = pos[h6];
						for ( i = 0; i<n1 ; i++ ) {
						tour[j] = h_tour[i]; 
						pos[h_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						for ( i = 0; i < n2 ; i++ ) {
						tour[j] = hh_tour[i];
						pos[hh_tour[i]] = j; 
						j++;
						if ( j >= n )
							j = 0;
						}
						tour[n] = tour[0];
					}
		      
					else if ( h == 2 ) {
			  
						/* copy part from pos[h1] to pos[h6]
						   (inverted) and from pos[h4] to pos[h5],
						   it remains the part from pos[h2] to
						   pos[h3] */
						// h1=c1  h6=s_c3
						j = pos[h1];
						h = pos[h6];
						i = 0;
						h_tour[i] = tour[j];
						n1 = 1;
						while ( j != h ) {
							/* 保存的是在tour中，逆序从h1(c1)->h6(s_c3)的点
								(但在h_tour中，是正序的)
							*/
							i++;
							j--;
							if ( j < 0 )
								j = n-1;
							h_tour[i] = tour[j];
							n1++;
						}

						j = pos[h4];
						h = pos[h5];
						i = 0;
						hh_tour[i] = tour[j];
						n2 = 1;
						while ( j != h ) {
							/* 保存的是在tour中，正序从h4(s_c2)->h5(c3)的点
								(但在hh_tour中，是正序的)
							*/
							i++;
							j++;
							if ( j >= n  )
								j = 0;
							hh_tour[i] = tour[j];
							n2++;
						}

						j = pos[h4];
						/* Now copy stored part from h_tour */
						for ( i = 0; i<n1 ; i++ ) {
							tour[j] = h_tour[i];
							pos[h_tour[i]] = j; 
							j++;
							if ( j >= n )
								j = 0;
						}
			  
						/* Now copy stored part from h_tour */
						for ( i = 0; i < n2 ; i++ ) {
							tour[j] = hh_tour[i];
							pos[hh_tour[i]] = j; 
							j++;
							if ( j >= n )
								j = 0;
						}
						tour[n] = tour[0];
					}
				}
				else if ( move_flag == 3 ) {
		      
				if ( pos_c3 < pos_c1 ) 
					n1 = pos_c1 - pos_c3;
				else
					n1 = n - (pos_c3 - pos_c1);
				if ( pos_c3 > pos_c2 ) 
					n2 = pos_c3 - pos_c2;
				else
					n2 = n - (pos_c2 - pos_c3);
				if ( pos_c2 > pos_c1 ) 
					n3 = pos_c2 - pos_c1;
				else
					n3 = n - (pos_c1 - pos_c2);
				/* n1: length h6 - h1, n2: length h4 - h5, n2: length h2 - h3 */
		      
				val[0] = n1; val[1] = n2; val[2] = n3; 
				/* Determine which is the longest part */
				h = 0;
				help = INT_MIN;
				for ( g = 0; g <= 2; g++) {
					if ( help < val[g] ) {
					help = val[g];
					h = g;
					}
				}
				/* order partial tours according length */
		      
				if ( h == 0 ) {
			  
					/* copy part from pos[h2] to pos[h3]
					   (inverted) and from pos[h4] to pos[h5]
					   it remains the part from pos[h6] to pos[h1] */
					j = pos[h3];
					h = pos[h2];
					i = 0;
					h_tour[i] = tour[j];
					n1 = 1;
					while ( j != h ) {
					i++;
					j--;
					if ( j < 0 )
						j = n-1;
					h_tour[i] = tour[j];
					n1++;
					}
			  
					j = pos[h2];
					h = pos[h5];
					i = pos[h4];
					tour[j] = h4;
					pos[h4] = j;
					while ( i != h ) {
					i++;
					if ( i >= n )
						i = 0;
					j++;
					if ( j >= n )
						j = 0;
					tour[j] = tour[i];
					pos[tour[i]] = j;
					}
					j++;
					if ( j >= n )
					j = 0;
					for ( i = 0; i < n1 ; i++ ) {
					tour[j] = h_tour[i];
					pos[h_tour[i]] = j; 
					j++;
					if ( j >= n )
						j = 0;
					}
					tour[n] = tour[0];
				}
				else if ( h == 1 ) {

					/* copy part from pos[h3] to pos[h2]
					   (inverted) and from  pos[h6] to pos[h1],
					   it remains the part from pos[h4] to pos[h5] */
					j = pos[h3];
					h = pos[h2];
					i = 0;
					h_tour[i] = tour[j];
					n1 = 1;
					while ( j != h ) {
					i++;
					j--;
					if ( j < 0  )
						j = n-1;
					h_tour[i] = tour[j];
					n1++;
					}

					j = pos[h6];
					h = pos[h1];
					i = 0;
					hh_tour[i] = tour[j];
					n2 = 1;
					while ( j != h ) {
					i++;
					j++;
					if ( j >= n )
						j = 0;
					hh_tour[i] = tour[j];
					n2++;
					}
			  
					j = pos[h6];
					for ( i = 0; i<n1 ; i++ ) {
					tour[j] = h_tour[i];
					pos[h_tour[i]] = j; 
					j++;
					if ( j >= n )
						j = 0;
					}

					for ( i = 0 ; i < n2 ; i++ ) {
					tour[j] = hh_tour[i];
					pos[hh_tour[i]] = j; 
					j++;
					if ( j >= n )
						j = 0;
					}
					tour[n] = tour[0];
				}
		      
				else if ( h == 2 ) {
			  
					/* copy part from pos[h4] to pos[h5]
					   (inverted) and from pos[h6] to pos[h1] (inverted)
					   it remains the part from pos[h2] to pos[h3] */
					j = pos[h5];
					h = pos[h4];
					i = 0;
					h_tour[i] = tour[j];
					n1 = 1;
					while ( j != h ) {
					i++;
					j--;
					if ( j < 0 )
						j = n-1;
					h_tour[i] = tour[j];
					n1++;
					}

					j = pos[h1];
					h = pos[h6];
					i = 0;
					hh_tour[i] = tour[j];
					n2 = 1;
					while ( j != h ) {
					i++;
					j--;
					if ( j < 0 )
						j = n-1;
					hh_tour[i] = tour[j];
					n2++;
					}

					j = pos[h4];
					/* Now copy stored part from h_tour */
					for ( i = 0; i< n1 ; i++ ) {
					tour[j] = h_tour[i];
					pos[h_tour[i]] = j; 
					j++;
					if ( j >= n )
						j = 0;
					}
					/* Now copy stored part from h_tour */
					for ( i = 0; i< n2 ; i++ ) {
					tour[j] = hh_tour[i];
					pos[hh_tour[i]] = j; 
					j++;
					if ( j >= n )
						j = 0;
					}
					tour[n] = tour[0];
				}
				}
				else {
				printf(" Some very strange error must have occurred !!!\n\n");
				exit(0);
				}
			}
			if (opt2_flag) {

				/* Now perform move */
				dlb[h1] = false; dlb[h2] = false;
				dlb[h3] = false; dlb[h4] = false;
				if ( pos[h3] < pos[h1] ) {
				help = h1; h1 = h3; h3 = help;
				help = h2; h2 = h4; h4 = help;
				}
				if ( pos[h3]-pos[h2] < n / 2 + 1) {
				/* reverse inner part from pos[h2] to pos[h3] */
				i = pos[h2]; j = pos[h3];
				while (i < j) {
					c1 = tour[i];
					c2 = tour[j];
					tour[i] = c2;
					tour[j] = c1;
					pos[c1] = j;
					pos[c2] = i;
					i++; j--;
				}
				}
				else {
				/* reverse outer part from pos[h4] to pos[h1] */
				i = pos[h1]; j = pos[h4];
				if ( j > i )
					help = n - (j - i) + 1;
				else 
					help = (i - j) + 1;
				help = help / 2;
				for ( h = 0 ; h < help ; h++ ) {
					c1 = tour[i];
					c2 = tour[j];
					tour[i] = c2;
					tour[j] = c1;
					pos[c1] = j;
					pos[c2] = i;
					i--; j++;
					if ( i < 0 )
					i = n - 1;
					if ( j >= n )
					j = 0;
				}
				tour[n] = tour[0];
				}
			}
	    }
	    else {
			dlb[c1] = true;
	    }
	} // end for ( l = 0 ; l < n ; l++ )
    }// end while(improveFlag)
	
}