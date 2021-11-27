// Author: Kwon Dohyun (Ajou Univ.)
// Last edited date: 11/27/2021

/*	
									DESCRIPTION

	This program calculates the structure of RNA's secondary structure that has maximum matching score and satisfy specific conditions.

	1. Watson-Crick pair and wobble pair is both allowed. If you want to set an arbitrary value, modify the argument parameters variable.
	2. Hairpin is allowed.
	3. Pseudoknot is allowed. But the section of the stem of the two Hairpin constituting the Pseudoknot is not allowed to overlap.
	4. Consider only when Hairpin and Pseudoknot are connected in series. This means that the structure surrounding the existing structure (nested structure) is not considered. (If you want this feature, you can modify it)

	TIME COMPLEXITY: O(N^5)
	SPACE COMPLEXITY: O(N^4)
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

#define MAX_SEQ_LENGTH 100

struct pk_info{
	int i,k,j;
};

// argument parameters
bool isWobblePairAllowed=false;
int wobblePairScore=1;
int watsonCrickPairScore=1;

// dp table
int lcs_dp[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int lcs_restore[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int pseudoknot[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
pk_info pseudoknot_restore[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int hairpin[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int hairpin_restore[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int structure[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];
int structure_restore[MAX_SEQ_LENGTH+1][MAX_SEQ_LENGTH+1];

vector<pii> restored;

// mischellaneous variables
int seq_len=-1;

void lcs(int s, int e, vector<vector<int>>& match){ // Time complexity - O(n^2)
	int len=e-s+1;
	for(int i=1;i<=len;i++){
		for(int j=1;j<=len;j++){
			if(i+j>len) break;
			lcs_dp[s][e][i][j]=lcs_dp[s][e][i-1][j];
			lcs_restore[s][e][i][j]=0;
			if(lcs_dp[s][e][i][j]<lcs_dp[s][e][i][j-1]){
				lcs_dp[s][e][i][j]=lcs_dp[s][e][i][j-1];
				lcs_restore[s][e][i][j]=1;
			}
			if(lcs_dp[s][e][i][j]<lcs_dp[s][e][i-1][j-1]+match[s+i-1][e-j+1]){
				lcs_dp[s][e][i][j]=lcs_dp[s][e][i-1][j-1]+match[s+i-1][e-j+1];
				lcs_restore[s][e][i][j]=2;
			}
		}
	}
}

void get_lcs_element(int s, int e, int l, int r){
	while(l&&r){
		if(lcs_restore[s][e][l][r]==0){ // (i-1, j) -> (i, j)
			l-=1;
		}
		else if(lcs_restore[s][e][l][r]==1){ // (i, j-1) -> (i, j)
			r-=1;
		}
		else{ // (i-1, j-1) + match -> (i, j)
			restored.push_back({s+l-1,e-r+1});
			l-=1; r-=1;
		}
	}
}

void restore(int l, int r){
	int len=r-l+1;
	//cout<<l<<" "<<r<<"\n";
	if(structure_restore[l][r]==-1){ // case of hairpin
		get_lcs_element(l,r,hairpin_restore[l][r],len-hairpin_restore[l][r]);
	}
	else if(structure_restore[l][r]==-2){ // case of pseudoknot
		auto [i,k,j]=pseudoknot_restore[l][r];
		get_lcs_element(l,j,i-l,j-k);
		get_lcs_element(i,r,k-i+1,r-j);
	}
	else{
		int k=structure_restore[l][r];
		if(l+1<=k&&l>=1){
			restore(l,k);	
		}
		if(k+2<=r&&k+1>=1){
			restore(k+1,r);
		}
	}
}

string get_dot_braket_notation(){
	string ret;
	ret.resize(seq_len, '.');
	
	sort(restored.begin(),restored.end());
	for(int i=0;i<(int)restored.size();i++){
		auto [l,r]=restored[i];
		l-=1;
		r-=1;		
		bool pk=false;
		for(int j=0;j<i;j++){
			auto [cl,cr]=restored[j];
			if(l<cr&&cr<r) {
				pk=true;
				break;
			}
		}
		if(!pk){
			ret[l]='(';
			ret[r]=')';
		}
		else{
			ret[l]='{';
			ret[r]='}';
		}
	}
	return ret;
}

int solve(string& seq){
	// make match matrix, if and only if (seq_i, seq_j) is Watson-Crick pair, set match[i][j] is true
	vector<vector<int>> match(seq_len+1,vector<int>(seq_len+1));
	for(int i=1;i<=seq_len;i++){
		for(int j=i+1;j<=seq_len;j++){
			if((seq[i]=='A'&&seq[j]=='U')||(seq[i]=='U'&&seq[j]=='A')||
				(seq[i]=='C'&&seq[j]=='G')||(seq[i]=='G'&&seq[j]=='C')){
				match[i][j]=match[j][i]=watsonCrickPairScore;
			}
			if(isWobblePairAllowed){ // if wobble pair?
				if((seq[i]=='G'&&seq[j]=='U')||(seq[i]=='U'&&seq[j]=='G')){
					match[i][j]=match[j][i]=wobblePairScore;
				}
			}
		}
	}

	// pre-compute lcs (longest common subsequence) - O(n^4)
	for(int i=1;i<=seq_len;i++){
		for(int j=i;j<=seq_len;j++){
			lcs(i,j,match);
		}
	}

	// get best hairpin (stem-loop) structure - O(n^3)
	for(int s=1;s<=seq_len;s++){
		for(int e=s+1;e<=seq_len;e++){
			int len=e-s+1;
			for(int k=0;k<=len;k++){
				if(hairpin[s][e]<lcs_dp[s][e][k][len-k]){
					hairpin[s][e]=lcs_dp[s][e][k][len-k];
					hairpin_restore[s][e]=k;	
				}				
			}
			//cout<<"s: "<<s<<" e: "<<e<<" dp: "<<hairpin[s][e]<<"\n";
		}
	}

	// get best pseudoknot structure - O(n^5)
	for(int s=1;s<=seq_len;s++){
		for(int e=s+3;e<=seq_len;e++){ // pseudoknot must have seq at least length 4
			int len=e-s+1;
			for(int i=s+1;i<e;i++){
				for(int j=s+1;j<e;j++){
					for(int k=i;k<j;k++){ // s - i - k - j - e
						if(pseudoknot[s][e]<lcs_dp[s][j][i-s][j-k]+lcs_dp[i][e][k-i+1][e-j]){
							pseudoknot[s][e]=lcs_dp[s][j][i-s][j-k]+lcs_dp[i][e][k-i+1][e-j];
							pseudoknot_restore[s][e]={i,k,j};
						}
					}
				}
			}
		}
	}

	// get best secondary structure
	for(int l=2;l<=seq_len;l++){
		for(int s=1;s+l-1<=seq_len;s++){
			int e=s+l-1;

			structure[s][e]=hairpin[s][e];
			structure_restore[s][e]=-1;

			if(structure[s][e]<pseudoknot[s][e]){
				structure[s][e]=pseudoknot[s][e];
				structure_restore[s][e]=-2;
			}

			for(int k=s;k<e;k++){
				if(structure[s][e]<structure[s][k]+structure[k+1][e]){
					structure[s][e]=structure[s][k]+structure[k+1][e];
					structure_restore[s][e]=k;
				}
			}
		}
	}
	return structure[1][seq_len];	
}

int main() {
	string seq; cin>>seq; // input RNA sequence with standard input
	seq_len=(int)seq.length();
	seq=" "+seq;

	int mfe=solve(seq); // calculate minimum free energy (MFE)

	cout<<" [ Input sequence ]\n";
	cout<<"RNA sequence:"<<seq<<"\n";
	cout<<"Length of sequence: "<<seq_len<<"\n";
	cout<<"\n";

	// print arguments
	cout<<" [ Arguments ]\n";
	cout<<"isWobblePairAllowed: ";
	if(isWobblePairAllowed) {
		cout<<"true\n";
		cout<<"Watson Crick Pair score: "<<watsonCrickPairScore<<"\n";
		cout<<"Wooble Pair score: "<<wobblePairScore<<"\n";
	}
	else cout<<"no\n";
	cout<<"\n";

	// get restored structure
	restore(1,seq_len);

	// print restored structure
	cout<<" [ Result ] \n";
	cout<<"Maximum matching score for mininum free energy (MFE): "<<mfe<<"\n";
	for(int i=1;i<(int)seq.size();i++) cout<<seq[i];
	cout<<"\n";
	cout<<get_dot_braket_notation()<<"\n";
	cout<<"\n";

	for(auto pi:restored){
		cout<<pi.first<<" "<<pi.second<<"\n";
	}
	return 0;
}

/*
hairpin example
ACGUGCCACGAUUCAACGUGGCACAG

pseudoknot example
UCGACUGUAAAGCGGCGACUUUCAGUCGCUCUUUUUGUCGCGCGC
*/