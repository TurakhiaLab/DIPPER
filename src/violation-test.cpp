#include <bits/stdc++.h>
using namespace std;
double eps=1e-5;
const int N=2000;
vector <pair<int,double>> v[N*2+5];
double td[N*2+5][N*2+5];
map <string,int> mp;
int cnt;
string s;
vector <int> oth;
void solve(int node,int l,int r){
    // cout<<node<<" "<<l<<" "<<r<<'\n';
    for(int j=l;j<=r;j++) cout<<s[j];puts("");
    for(int i=l;i<=r;i++){
        // cout<<i<<'\n';
        if(s[i]=='('){
            int id=++cnt;
            solve(id,i+1,oth[i]-1);
            i=oth[i]+2;
            string temp;
            while(i<=r&&((s[i]>='0'&&s[i]<='9')||s[i]=='.'||s[i]=='e'||s[i]=='-')) temp+=s[i++];
            v[node].push_back({id,stod(temp)});
            v[id].push_back({node,stod(temp)});
            // cout<<node<<" "<<id<<"\n";
            assert(s[i]==','||s[i]==')');
        }
        else{
            int id=++cnt;
            string temp;
            while(s[i]!=':') temp+=s[i++];
            mp[temp]=id;
            string dig;
            i++;
            while(i<=r&&((s[i]>='0'&&s[i]<='9')||s[i]=='.'||s[i]=='e'||s[i]=='-')) dig+=s[i++];
            v[node].push_back({id,stod(dig)});
            v[id].push_back({node,stod(dig)});
            // cout<<node<<" "<<id<<'\n';
        }
    }
}
string name[N+5];
double dis[N+5][N+5];
int main(){
    string file1 = "/home/zec022@AD.UCSD.EDU/phylo-accel/test/phastsim_dataset/estimated_tree.nwk";
    string file2 = "/data/zec022/distance_matrix.phy";
    ifstream inputFile1(file1);
    getline(inputFile1,s);
    oth.resize(s.size());
    vector <int> stk;
    for(int i=0;i<s.size();i++){
        if(s[i]=='(' ) stk.push_back(i);
        else if(s[i]==')') oth[stk.back()]=i,stk.pop_back();
    }
    solve(0,1,int(s.size())-2);
    for(int i=0;i<=cnt;i++){
        for(int j=0;j<=cnt;j++) td[i][j]=-1;
        td[i][i]=0;
        queue <int> q;
        q.push(i);
        while(!q.empty()){
            int temp=q.front();
            q.pop();
            for(auto t:v[temp])
                if(td[i][t.first]==-1)
                    td[i][t.first]=td[i][temp]+t.second,q.push(t.first);
        }
    }
    ifstream inputFile2(file2);
    int n;
    inputFile2>>n;
    for(int i=0;i<n;i++){
        inputFile2>>name[i];
        for(int j=0;j<n;j++) inputFile2>>dis[i][j];
    }
    for(int i=0;i<n;i++) cout<<name[i]<<' '<<mp[name[i]]<<'\n';
    int tot=n*n,bad=0;
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            if(td[mp[name[i]]][mp[name[j]]]<dis[i][j]+eps)
                bad++;
            // cout<<dis[i][j]<<" "<<td[mp[name[i]]][mp[name[j]]]<<'\n';
    cout<<double(bad)/tot<<'\n';
}