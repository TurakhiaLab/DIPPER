#include <bits/stdc++.h>
using namespace std;
int main(){
    int n;
    cin>>n;
    vector <string> name;
    name.resize(n);
    vector <vector<double>> dist;
    dist.resize(n);
    for(int i=0;i<n;i++){
        cin>>name[i];
        dist[i].resize(n);
        for(int j=0;j<n;j++) scanf("%lf",&dist[i][j]);
    }
    int bd=n/10;
    for(int i=0;i<n;i++){
        cout<<'\t'<<name[i];
    }
    cout<<'\n';
    for(int i=bd;i<n;i++){
        cout<<name[i]<<'\t';
        for(int j=0;j<=n;j++) printf("%lf\t",dist[i][j]);
        cout<<'\n';
    }
}