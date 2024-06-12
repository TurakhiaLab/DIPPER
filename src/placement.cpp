/*
Input format: phylip format distance matrix, with two additional parameters on the first line

Sample:
N bd T
Name[0] dis[0][0] dis[0][1] ... dis[0][N-1]
...
Name[N-1] dis[N-1][0] dis[N-1][1] ... dis[N-1][N-1]

bd should between 2 to N, preferably 3, indicate we will build a NJ tree of size bd at the beginning
T could be any positive value, indicating that the program will generate phylogenetic trees for T times
All T trees being generated are based on the given distance matrix,
but with different placement orders (this would help testing the average value)

Output format:
N trees, each of them being quoted, and separated by space

Sample:
'Tree0' 'Tree1' 'Tree2'

 */
#include <bits/stdc++.h>
using namespace std;
double dis[5005][5005],odis[5005][5005];
string oname[5005],name[5005];
double u[5005];
int f[5005],id[5005],rid[5005];
vector <pair<int,double>> v[5005];
string ttt;
void print(int x,int from) {
    //cout<<"#"<<x<<endl;
    if (v[x].size()>1) {
        int cnt=0;
//        printf("(");
        ttt+='(';
        for (int i = 0; i < v[x].size(); i++) {
            if(v[x][i].first==from) continue;
            cnt++;
            print(v[x][i].first,x);
//            printf(":%.5lf%c", v[x][i].second, cnt==2 ? ')' : ',');
            ttt+=':';
            ttt+=to_string(v[x][i].second);
            ttt+=cnt==2?')':',';
        }
    }
    else{
//        cout<<name[x];
        ttt+=name[x];
    }
}
double mnerror;
pair <int,int> bestpos;
pair <double,double> bl;
double tar_bl;
int under[5005];
int flag[5005];
void color(int x,int from){
    under[x]=1;
    for(auto t:v[x]) if(t.first!=from) color(t.first,x);
}
void dfs(int x,int from,int tot){
    for(auto t:v[x]) {
        if(t.first==from) continue;
        dfs(t.first,x,tot);
        if(t.second==0) continue;
        for(int i=0;i<tot;i++) under[i]=0;
        color(t.first,x);
        double disunder=0,disup=0;
        queue <pair<int,double>> q;
        int vis[505]={},td[505]={};
        double tdis[505]={};
        q.push({x,0}),td[x]=1,tdis[x]=0;
        vector <int> idup,idunder;
        int tc=0;
        while(!q.empty()&&tc<=5){
            auto temp=q.front();
            q.pop();
            if(temp.first<tot) disup=max(disup,dis[temp.first][tot]-temp.second),tc++,idup.push_back(temp.first);
            vis[temp.first]=1;
            for(auto nxt:v[temp.first])
                if(!vis[nxt.first]&&nxt.first!=t.first)
                    q.push({nxt.first,temp.second+nxt.second}),td[nxt.first]=td[temp.first]+1,tdis[nxt.first]=temp.second+nxt.second;
        }
        while(!q.empty()) q.pop();
        tc=0;
        q.push({t.first,0}),td[t.first]=1;
        while(!q.empty()&&tc<=5){
            auto temp=q.front();
            q.pop();
            if(temp.first<tot) disunder=max(disunder,dis[temp.first][tot]-temp.second),tc++,idunder.push_back(temp.first);
            vis[temp.first]=1;
            for(auto nxt:v[temp.first])
                if(!vis[nxt.first]&&nxt.first!=x)
                    q.push({nxt.first,temp.second+nxt.second}),td[nxt.first]=td[temp.first]+1,tdis[nxt.first]=temp.second+nxt.second;
        }
        double additional_dis=max(0.0,(disup+disunder-t.second)/2);
        disunder-=additional_dis,disup-=additional_dis;
        if(disunder<0) disunder=0;
        if(disup<0) disup=0;
        //if(aveup+aveunder==0) continue;
        double dis_to_under=disunder,dis_to_up=disup;
        double error=additional_dis;
//        for(auto id:idup) error+=(additional_dis+dis_to_up+tdis[id])/(1<<td[id]);
//        for(auto id:idunder) error+=(additional_dis+dis_to_under+tdis[id])/(1<<td[id]);
        //accum(t.first,x,dis_to_under+additional_dis,error,tot);
        //accum(x,t.first,dis_to_up+additional_dis,error,tot);
//        cout<<x<<" "<<t.first<<" "<<t.second<<endl;
//        cout<<dis_to_up<<" "<<additional_dis<<endl;
//        cout<<error<<endl<<endl;
        if(error<mnerror){
            mnerror=error;
            bestpos={x,t.first};
            bl={dis_to_up,additional_dis};
            tar_bl=t.second;
        }
    }
}
double dep[5005];
void find_dep(int x,int from){
    for(auto t:v[x])
        if(t.first!=from)
            dep[t.first]=dep[x]+t.second,find_dep(t.first,x);
}
vector<pair<double,string>> output;
void get_id(int x,int from,vector<int>&ids){
    ids.push_back(x);
    for(auto t:v[x])
        if(t.first!=from)
            get_id(t.first,x,ids);
}
map <vector<int>,double> qryans;
int n;
double query(int x,int fx,int y,int fy){
    vector <int> cx,cy;
    get_id(x,fx,cx),get_id(y,fy,cy);
    double ans=0;
    for(auto idx:cx)
        for(auto idy:cy)
            ans+=dis[idx][idy];
    return ans/cx.size()/cy.size();


//    vector <int> tp;
//    tp.push_back(x),tp.push_back(fx),tp.push_back(y),tp.push_back(fy);
//    if(qryans.count(tp)) return qryans[tp];
//    double ans=0;
//    if(x<n){
//        if(y<n) ans=dis[x][y];
//        else for(auto t:v[y]) if(t.first!=fy) ans+=query(x,fx,t.first,y)/2;
//    }
//    else for(auto t:v[x]) if(t.first!=fx) ans+=query(t.first,x,y,fy)/2;
//    return qryans[tp]=ans;
}
void assign_weight(int x,int from){
    for(auto &t:v[x])
        if(t.first!=from){
            assign_weight(t.first,x);
            if(t.first<n){
                t.second=0;
                vector <int> spid;
                int xx,oth;
                if(v[x].size()==2){
                    for(auto tt:v[x]) if(tt.first!=t.first) xx=tt.first;
                    oth=x;
                }
                else xx=x,oth=t.first;
                for(auto tt:v[xx]){
                    if(tt.first!=oth) spid.push_back(tt.first);
                }
                assert(spid.size()==2);
                t.second+=query(t.first,x,spid[0],xx);
                t.second+=query(t.first,x,spid[1],xx);
                t.second-=query(spid[0],xx,spid[1],xx);
                t.second/=2;
                if(xx!=x) t.second/=2;
            }
            else{
                vector <int> low,high;
                for(auto tt:v[t.first]) if(tt.first!=x) low.push_back(tt.first);
                int xx,oth;
                if(v[x].size()==2){
                    for(auto tt:v[x]) if(tt.first!=t.first) xx=tt.first;
                    oth=x;
                }
                else xx=x,oth=t.first;
                if(v[xx].size()==1){
                    t.second=0;
                    vector <int> spid;
                    int sub;
                    for(auto tt:v[x]) if(tt.first!=t.first) sub=tt.first;
                    int xx,oth;
                    if(v[x].size()==2){
                        for(auto tt:v[x]) if(tt.first!=sub) xx=tt.first;
                        oth=x;
                    }
                    else xx=x,oth=sub,assert(0);
                    for(auto tt:v[xx]){
                        if(tt.first!=oth) spid.push_back(tt.first);
                    }
                    assert(spid.size()==2);
                    t.second+=query(sub,x,spid[0],xx);
                    t.second+=query(sub,x,spid[1],xx);
                    t.second-=query(spid[0],xx,spid[1],xx);
                    t.second/=2;
                    if(xx!=x) t.second/=2;
                    continue;
                }
                for(auto tt:v[xx]) if(tt.first!=oth) high.push_back(tt.first);
                assert(low.size()==2);
                if(high.size()!=2){
                    cerr<<x<<endl;
                    cerr<<xx<<" "<<oth<<'\n';
                    cerr<<v[xx].size()<<' '<<high.size()<<'\n';
                }
                assert(high.size()==2);
                vector <int> va,vb,vc,vd;
                get_id(low[0],t.first,va),get_id(low[1],t.first,vb);
                get_id(high[0],xx,vc),get_id(high[1],xx,vd);
                double lam;
                lam=double(va.size()*vd.size()+vb.size()*vc.size())/(va.size()+vb.size())/(vc.size()+vd.size());
//                lam=0.5;
                t.second=0;
                t.second+=lam*(query(low[0],t.first,high[0],xx)+query(low[1],t.first,high[1],xx));
                t.second+=(1-lam)*(query(low[0],t.first,high[1],xx)+query(low[1],t.first,high[0],xx));
                t.second-=query(low[0],t.first,low[1],t.first)+query(high[0],x,high[1],xx);
                t.second/=2;
                if(xx!=x) t.second/=2;
            }
        }
}

int main() {
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);
    int bd,totT;
    cin>>n>>bd>>totT;
    for(int i=0;i<n;i++) {
        cin>>oname[i];
        for (int j = 0; j < n; j++)
            scanf("%lf", &odis[i][j]);
    }
    mt19937 rnd(time(NULL));
    for(int T=0;T<totT;T++) {
        cerr<<T<<'\n';
        for (int i = 0; i < n; i++) id[i] = i, rid[i] = i,f[i]=0,v[i*3].clear(),v[i*3+1].clear(),v[i*3+2].clear(),u[i]=0,dep[i]=0,flag[i]=under[i]=0;
        shuffle(rid, rid + n, rnd);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                dis[i][j] = odis[rid[i]][rid[j]];
        for(int i=0;i<n;i++) name[i]=oname[rid[i]];
        int cnt = n;
        for (int i = 0; i < bd - 1; i++) {
            if (i == bd - 2) {
                vector<int> vv;
                for (int j = 0; j < bd; j++)
                    if (!f[j]) vv.push_back(j);
                double dd = dis[vv[0]][vv[1]];
                v[cnt].push_back({id[vv[0]], dd / 2}), v[cnt].push_back({id[vv[1]], dd / 2});
                v[id[vv[0]]].push_back({cnt,dd/2}),v[id[vv[1]]].push_back({cnt,dd/2});
                cnt++;
                continue;
            }
            for (int j = 0; j < bd; j++) {
                u[j] = 0;
                for (int k = 0; k < bd; k++)
                    if (f[j] == 0 && f[k] == 0 && j != k)
                        u[j] += dis[j][k];
                u[j] /= (bd - i - 2);
            }
            double mn = 10000;
            int x = -1, y = -1;
            for (int j = 0; j < bd; j++) {
                if (f[j]) continue;
                for (int k = 0; k < bd; k++) {
                    if (j == k || f[k]) continue;
                    double val = dis[j][k] - u[j] - u[k];
                    if (val < mn) mn = val, x = j, y = k;
                }
            }
            double blx = (dis[x][y] + u[x] - u[y]) * 0.5;
            double bly = dis[x][y] - blx;
            if (blx < 0) bly += blx, blx = 0;
            if (bly < 0) blx += bly, bly = 0;
            v[cnt].push_back({id[x], blx}), v[cnt].push_back({id[y], bly});
            v[id[x]].push_back({cnt, blx}), v[id[y]].push_back({cnt, bly});
            //cerr << i<<" "<<id[x] << " " << id[y] << " ";
            //printf("%.10lf %.10lf\n",mn, dis[x][y]);
            f[y] = 1, id[x] = cnt++;
            for (int j = 0; j < bd; j++)
                if (j != x && j != y && !f[j]) {
                    double val = (dis[x][j] + dis[y][j] - dis[x][y]) / 2.0;
                    //if(val<0) val=0;
                    dis[x][j] = dis[j][x] = val;
                }
        }
        int rt = cnt - 1;
        for (int i = bd; i < n; i++) {
//            cerr<<T<<" "<<i<<'\n';
            for (int j = 0; j < i; j++) flag[j]=1;
            mnerror = 1e9;
            dfs(rt, -1, i);
//        cout<<mnerror<<endl;
//        cout<<bestpos.first<<" "<<bestpos.second<<endl;
//        cout<<bl.first<<" "<<bl.second<<endl;
            v[bestpos.first].erase(
                    std::remove_if(v[bestpos.first].begin(), v[bestpos.first].end(), [&](const std::pair<int, int> &p) {
                        return p.first == bestpos.second;
                    }), v[bestpos.first].end());
            v[bestpos.second].erase(std::remove_if(v[bestpos.second].begin(), v[bestpos.second].end(),
                                                   [&](const std::pair<int, int> &p) {
                                                       return p.first == bestpos.first;
                                                   }), v[bestpos.second].end());
            v[bestpos.first].push_back({cnt, bl.first}), v[bestpos.second].push_back({cnt, tar_bl - bl.first});
            v[cnt].push_back({bestpos.first, bl.first}), v[cnt].push_back({bestpos.second, tar_bl - bl.first});
            v[cnt].push_back({i, bl.second}), v[i].push_back({cnt, bl.second});
            cnt++;
        }
        assign_weight(rt,-1);
        ttt="";
        ttt+='\'';
//        printf("\'");
        print(rt, -1);
        ttt+=";\' ";
//    for(int i=bd;i<n;i++){
//        dfs(rt);
//    }
        for (int i = 0; i < n; i++) flag[i] = 1;
        double sum = 0;
        for (int i = 0; i < n; i++) for (auto t:v[i]) sum += t.second;
//        cout << sum << '\n';
        output.push_back({sum,ttt});
    }
    sort(output.begin(),output.end());
//    for(int i=0;i<10;i++) cout<<output[i].second;
    for(int i=0;i<totT;i++){
        double mn=1e9;
        string ans;
        for(int j=i;j<i+1;j++) if(output[j].first<mn){
                mn=output[j].first,ans=output[j].second;
            }
        cout<<ans;
    }
    return 0;
}
