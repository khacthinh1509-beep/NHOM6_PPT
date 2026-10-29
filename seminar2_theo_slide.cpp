#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

using So = long double;
using Day = vector<So>;
using Bang = vector<Day>;
// He so da thuc luu theo bac TANG: p[0] + p[1]*x + ...
Day nhanNhanTu(const Day& p, So a) {
    Day q(p.size()+1,0);
    for (size_t i=0;i<p.size();++i) {
        q[i]-=a*p[i]; q[i+1]+=p[i];
    }
    return q; // p(x)*(x-a)
}
So giaTri(const Day& p, So x) {
    So y=0;
    for (int i=static_cast<int>(p.size())-1;i>=0;--i) y=y*x+p[i];
    return y;
}
// Slide chuong 2: p_n(x) = a_0*x^n + ... + a_n.
// Mang noi bo luu bac tang; khi in dao thu tu de dung voi slide.
void inDaThuc(const char* ten, const Day& heSo) {
    int n=static_cast<int>(heSo.size())-1;
    cout << ten << "(x) = ";
    bool dau=true;
    for(int k=n;k>=0;--k) {
        So v=heSo[k];
        if(abs(v)<1e-12L) continue; // Chi an sai so lam tron khi in.
        if(!dau) cout << (v<0 ? " - " : " + ");
        else if(v<0) cout << "-";
        cout << abs(v);
        if(k>0) cout << "*x";
        if(k>1) cout << "^" << k;
        dau=false;
    }
    if(dau) cout << "0";
    cout << '\n';
    cout << "He so theo bac giam (a_0,...,a_" << n << "): ";
    for(int k=n;k>=0;--k)
        cout << (abs(heSo[k])<1e-12L ? 0.0L : heSo[k]) << "  ";
    cout << '\n';
}
void doiChieu(const char* ten, So may, So tay) {
    cout << ten << ": may=" << may << "; tay=" << tay
              << "; lech=" << scientific << abs(may-tay)
              << fixed << '\n';
    if (abs(may-tay)>1e-8L) throw runtime_error("Doi chieu khong khop");
}
// n la bac toi da; co n+1 moc x_0,...,x_n nhu slide.
Day lagrange(const Day& x, const Day& y) {
    int n=static_cast<int>(x.size())-1;
    Day p_n(n+1,0);
    for(int i=0;i<=n;++i) {
        Day l_i{1}; So mau=1;
        for(int j=0;j<=n;++j) if(j!=i) {
            l_i=nhanNhanTu(l_i,x[j]);
            mau*=x[i]-x[j];
        }
        if(mau==0) throw runtime_error("Cac moc phai khac nhau");
        for(int k=0;k<=n;++k) p_n[k]+=y[i]*l_i[k]/mau;
    }
    return p_n;
}
// Bai 4 co moc KHONG cach deu: dung ty sai phan f[x_i,...,x_(i+k)].
// Khong duoc thay bang Delta^k y_i cua cong thuc moc cach deu.
// Khi moc cach deu: f[x_i,...,x_(i+k)] = Delta^k y_i/(k!*h^k).
Day newton(const Day& x,const Day& y) {
    int n=static_cast<int>(x.size())-1;
    Bang T(n+1,Day(n+1,0));
    for(int i=0;i<=n;++i) T[i][0]=y[i];
    for(int k=1;k<=n;++k)
        for(int i=0;i<=n-k;++i)
            T[i][k]=(T[i+1][k-1]-T[i][k-1])/(x[i+k]-x[i]);
    cout << "Moc khong cach deu; T[i][k]=f[x_i,...,x_(i+k)]\n";
    cout << "Bang ty sai phan (cac cot tu cap 0 den cap n):\n";
    for(int i=0;i<=n;++i) {
        for(int k=0;k<=n-i;++k) cout << T[i][k] << "  ";
        cout << '\n';
    }
    Day L_n{T[0][n]};
    for(int k=n-1;k>=0;--k) {
        L_n=nhanNhanTu(L_n,x[k]); L_n[0]+=T[0][k];
    }
    return L_n;
}
void bai1() {
    cout << "\nBAI 1 - LAGRANGE\n";
    Day p=lagrange({0,1,3,5},{1,2,1,4});inDaThuc("p_3",p);
    doiChieu("p_3(2)",giaTri(p,2),8.0L/5);
}
void bai2() {
    cout << "\nBAI 2 - LAGRANGE VA PHAN THUC\n";
    Day x{-4,-3,-1,0},y;for(So t:x)y.push_back(t*t*t+t*t-10);
    inDaThuc("p_3",lagrange(x,y));Day r;
    for(int i=0;i<4;++i){So m=1;for(int j=0;j<4;++j)if(i!=j)m*=x[i]-x[j];r.push_back(y[i]/m);}
    doiChieu("A",r[3],-5.0L/6);doiChieu("B",r[2],5.0L/3);
    doiChieu("C",r[1],-14.0L/3);doiChieu("D",r[0],29.0L/6);
}
void bai3() {
    cout << "\nBAI 3 - TONG LAP PHUONG\n";
    Day p=lagrange({0,1,2,3,4},{0,1,9,36,100});inDaThuc("p_4",p);
    cout << "S_n=p_4(n)=[n*(n+1)/2]^2\n";
    for(int n:{1,5,10,20}) {
        So s=0;for(int k=1;k<=n;++k)s+=So(k)*k*k;
        cout << "n=" << n << "  ";doiChieu("S(n)",giaTri(p,n),s);
    }
}
void bai4() {
    cout << "\nBAI 4 - NOI SUY NEWTON\n";
    Day p=newton({-1,0,3,6,7},{3,-6,39,822,1611});inDaThuc("L_4",p);
    doiChieu("L_4(0.25)",giaTri(p,.25L),-1467.0L/256);
}
int main() {
    cout << fixed << setprecision(10);
    cout << "BAI 9 SEMINAR 2 - CHUONG 2\n";
    cout << "Cac bai 1, 2, 3, 4: noi suy va sai so noi suy.\n";
    cout << "Ky hieu: p_n, l_i (Lagrange); L_n (Newton).\n";
    try {
        bai1(); bai2(); bai3(); bai4();
    } catch(const exception& e) {
        cerr << "Loi: " << e.what() << '\n';
        return 1;
    }
    cout << "\nHoan thanh 4 bai thuoc Chuong 2.\n";
    return 0;
}
