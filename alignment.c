#include<stdio.h>
#include<stdbool.h>
#include<stdlib.h>
#include<time.h> 
#define M 15 
#define N 15 
#define SEED_1 200
#define SEED_2 6

typedef struct{
	int i;
	int j;
	int opt;
	struct NODE* pre;
}NODE;

void print_sequence(char* str, int length);
bool get_sequence(char* str, int length, int seed);
int min(int t1, int t2, int t3);
int opt_cost_recursion(char* x, char* y, int i, int j, int m, int n);
int opt_cost_dynamic(char* x, char* y, int m, int n, NODE** node);
void valuate_recursion(char* x, char* y, int m, int n);
void valuate_dynamic(char* x, char* y, int m, int n);

int main() {
	int opt;
	clock_t start, end;
	double dur;
	int m = M;
	int n = N;
	char* x = (char*)malloc((m+1)*sizeof(char));
	char* y = (char*)malloc((n+1)*sizeof(char));
	if(get_sequence(x, m, SEED_1) && get_sequence(y, n, SEED_2)) {
		printf("DNA sequences: \n");
		print_sequence(x, m);
		print_sequence(y, n);
		printf("\n-----Recursion result-----\n");
		valuate_recursion(x, y, m, n);
		printf("-----Dynamic planning result-----\n");
		valuate_dynamic(x, y, m, n);
	}
	free(x);
	free(y);
	return 0;
}

void print_sequence(char* str, int length) {
	int i;
	for(i = 0; i < length; i++) {
		printf("%c", str[i]);
	}
	printf("\n");
	return;
}

bool get_sequence(char* str, int length, int seed) {
	int i;
	int randnum;
	srand((unsigned int)seed);
	for(i = 0; i < length; i++) {
		randnum = rand()%4;
		if(randnum == 0) str[i] = 'A';
		else if(randnum == 1) str[i] = 'T';
		else if(randnum == 2) str[i] = 'C';
		else if(randnum == 3) str[i] = 'G';
		else return false;
	}
	str[i] = '\0';
	return true;
}

int min(int t1, int t2, int t3) {
	int min_num;
	min_num = (t1 < t2) ? t1 : t2;
	min_num = (t3 < min_num) ? t3 : min_num;
	return min_num;
}

int opt_cost_recursion(char* x, char* y, int i, int j, int m, int n) {
	int penalty;
	if(i == m) {
		return 2*(n-j);
	}
	else if(j == n) {
		return 2*(m-i);
	}
	else {
		penalty = (x[i] == y[j]) ? 0 : 1;
		return min(opt_cost_recursion(x, y, i+1, j+1, m, n) + penalty, 
					opt_cost_recursion(x, y, i, j+1, m, n) + 2,
					opt_cost_recursion(x, y, i+1, j, m, n) + 2);
	}
}


int opt_cost_dynamic(char* x, char* y, int m, int n, NODE** node) {
	int i, j, k;
	int penalty;
	int opt_temp;
	NODE* node_temp;

	// initialize termination node[0..m][n], node[m][0..n]
	for(i = m; i >= 0; i--) {
		node[i][n].opt = 2 * (m - i);
		node[i][n].pre = NULL;
	}
	for(j = n; j >= 0; j--) {
		node[m][j].opt = 2 * (n - j);
		node[m][j].pre = NULL;
	}
	// fill node[0..(m-1)][0..(n-1)]
	for(i = m-1; i >= 0; i--) {
		for(j = n-1; j>= 0; j--) {
			penalty = (x[i] == y[j]) ? 0 : 1;
			if(node[i][j+1].opt+2 < node[i+1][j].opt+2) {
				opt_temp = node[i][j+1].opt + 2;
				node_temp = &(node[i][j+1]);
			} else {
				opt_temp = node[i+1][j].opt + 2;
				node_temp = &(node[i+1][j]);
			}
			if(node[i+1][j+1].opt+penalty < opt_temp) {
				opt_temp = node[i+1][j+1].opt + penalty;
				node_temp = &(node[i+1][j+1]);
			}
			node[i][j].opt = opt_temp;
			node[i][j].pre = (struct NODE*)node_temp;
		}
	}
	return node[0][0].opt;
}

void valuate_recursion(char* x, char* y, int m, int n) {
	clock_t start, end;
	double dur;
	int opt;
	start = clock();
	opt = opt_cost_recursion(x, y, 0, 0, m, n);
	end = clock();
	dur = (double)(end - start);
	printf("Optimal Cost: %d\n", opt);
	printf("Use Time: %f s\n\n", (dur/CLOCKS_PER_SEC));
	return;
}

void valuate_dynamic(char* x, char* y, int m, int n) {
	clock_t start, end;
	double dur;
	int i, j, k, opt;
	int p_i, p_j;
	int aligned_length;
	char* aligned_x;
	char* aligned_y;
	NODE *temp;
	NODE **node;
	
	start = clock();
	// allocate space NODE node[m+1][n+1]
	node = (NODE**)malloc((m+1)*sizeof(NODE*));
	for(k = 0; k <= m; k++) {
		node[k] = (NODE*)malloc((n+1)*sizeof(NODE));
	}
	// initialize all nodes
	for(i = 0; i <= m; i++) {
		for(j = 0; j <= n; j++) {
			node[i][j].i = i;
			node[i][j].j = j;
			node[i][j].opt = 0;
			node[i][j].pre = NULL;
		}
	}
	opt = opt_cost_dynamic(x, y, m, n, node);
	end = clock();
	dur = (double)(end - start);

	// allocate the optimal alignment sequences' space
	aligned_length = m + n;
	aligned_x = (char*)malloc((aligned_length+1)*sizeof(char));
	aligned_y = (char*)malloc((aligned_length+1)*sizeof(char));
	// initialize
	for(k = 0; k < aligned_length; k++) {
		aligned_x[k] = '\0';
		aligned_y[k] = '\0';
	}
	// fill the optimal alignment sequences
	temp = &node[0][0];
	k = 0;
	while(temp->pre != NULL) {
		i = temp->i;
		j = temp->j;
		temp = (NODE*)temp->pre;
		p_i = temp->i;
		p_j = temp->j;
		if((p_i == (i+1)) && (p_j == (j+1))) {
			aligned_x[k] = x[i];
			aligned_y[k] = y[j];
		} else if((p_i == i) && (p_j == (j+1))) {
			aligned_x[k] = '-';
			aligned_y[k] = y[j];
		} else if((p_i == (i+1)) && (p_j == j)) {
			aligned_x[k] = x[i];
			aligned_y[k] = '-';
		} else {
			printf("Alignment Error!");
			exit(1);
		}
		k++;
	}
	if(p_i < m && p_j < n) { // pre chain break incorrectly
		printf("Alignment Error!");
		exit(2);
	} else if(p_i < m) {
		for(i = p_i; i < m; i++) {
			aligned_x[k] = x[i];
			aligned_y[k] = '-';
			k++;
		}
	} else if(p_j < n) {
		for(j = p_j; j < n; j++) {
			aligned_x[k] = '-';
			aligned_y[k] = y[j];
			k++;
		}
	}
	printf("Optimal Cost: %d\n", opt);
	printf("Optimal Alignment: \n");
	print_sequence(aligned_x, aligned_length);
	print_sequence(aligned_y, aligned_length);
	printf("Use Time: %f s\n\n", (dur/CLOCKS_PER_SEC));
	// free
	free(aligned_x);
	free(aligned_y);
	for(k = 0; k <= m; k++) {
		free(node[k]);
	}
	free(node);
	return;
}
