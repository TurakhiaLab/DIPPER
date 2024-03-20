#include <bits/stdc++.h>

using namespace std;
int main(){
    mt19937 myrnd(time(NULL));
    int n;
    cin>>n;
    double **dis = new double*[n];
    for (int i = 0; i < n; ++i) {
        dis[i] = new double[n];
    }
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            dis[i][j]=dis[j][i]=double(myrnd()%100000)/100000;
    cout<<n<<"\n";
    for(int i=0;i<n;i++){
        printf("Node_%d\t",i);
        for(int j=0;j<n;j++){
            printf("%.6lf%c",dis[i][j],j+1==n?'\n':'\t');
        }
    }
    return 0;
}