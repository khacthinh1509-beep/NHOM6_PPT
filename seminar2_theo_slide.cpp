#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
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
void inMoHinh(const Day& heSo) {
    So a=heSo[0], b=heSo[1];
    if(heSo.size()==2)
        cout << "y=a+b*x: a=" << a << "; b=" << b << '\n';
    else {
        So c=heSo[2];
        cout << "y=a+b*x+c*x^2: a=" << a
                  << "; b=" << b << "; c=" << c << '\n';
    }
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
So lagrangeTai(const Day& x,const Day& y,So t) {
    So p_n=0;
    int n=static_cast<int>(x.size())-1;
    for(int i=0;i<=n;++i) {
        So l_i=1;
        for(int j=0;j<=n;++j)
            if(i!=j) l_i*=(t-x[j])/(x[i]-x[j]);
        p_n+=y[i]*l_i;
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
Day gauss(Bang A,Day b) {
    int n=static_cast<int>(b.size());
    for(int k=0;k<n;++k) {
        int r=k;
        for(int i=k+1;i<n;++i) if(abs(A[i][k])>abs(A[r][k])) r=i;
        if(abs(A[r][k])<1e-18L) throw runtime_error("He suy bien");
        swap(A[k],A[r]);swap(b[k],b[r]);
        for(int i=k+1;i<n;++i) {
            So m=A[i][k]/A[k][k];
            for(int j=k;j<n;++j) A[i][j]-=m*A[k][j];
            b[i]-=m*b[k];
        }
    }
    Day c(n,0);
    for(int i=n-1;i>=0;--i) {
        So s=b[i];for(int j=i+1;j<n;++j)s-=A[i][j]*c[j];
        c[i]=s/A[i][i];
    }
    return c;
}
Day binhPhuongToiThieu(const Day& x,const Day& y,int bac) {
    Bang A(bac+1,Day(bac+1,0));Day b(bac+1,0);
    // A[k][j] = tong x_i^(k+j), b[k] = tong x_i^k*y_i.
    for(size_t i=0;i<x.size();++i) {
        Day mu(2*bac+1,1);
        for(int j=1;j<=2*bac;++j)mu[j]=mu[j-1]*x[i];
        for(int k=0;k<=bac;++k) {
            b[k]+=mu[k]*y[i];
            for(int j=0;j<=bac;++j)A[k][j]+=mu[k+j];
        }
    }
    return gauss(A,b);
}
// Chuong 3: epsilon_i = y_i - f(x_i); S = tong epsilon_i^2.
// RMSE la chi so bo sung. Chi so mang C++ bat dau tu 0.
template<class Ham>
So inSaiSo(const char* ten,const Day& x,const Day& y,Ham f) {
    So S=0;
    for(size_t i=0;i<x.size();++i) {
        So epsilon_i=y[i]-f(x[i]);
        S+=epsilon_i*epsilon_i;
    }
    cout << ten << ": S=" << S
              << "; RMSE=" << sqrt(S/x.size()) << '\n';
    return S;
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
void bai5() {
    cout << "\nBAI 5 - BON DANG XAP XI\n";
    Day x{0,1,2,3,4,5,6,7},y{1.4L,1.3L,1.4L,1.1L,1.3L,1.8L,1.6L,2.3L};
    for(int d:{1,2}) {
        Day p=binhPhuongToiThieu(x,y,d);inMoHinh(p);
        inSaiSo(d==1?"Bac 1":"Bac 2",x,y,[&](So t){return giaTri(p,t);});
    }
    Day Y;for(So t:y)Y.push_back(log(t));
    Day heSo=binhPhuongToiThieu(x,Y,1);
    So A=heSo[0], B=heSo[1];
    So a=exp(A), b=B;
    cout << "Y=ln(y), X=x; Y=A+B*X; A=" << A << "; B=" << B << '\n';
    cout << "y=a*exp(b*x): a=" << a << "; b=" << b << '\n';
    inSaiSo("Ham mu",x,y,[&](So t){return a*exp(b*t);});
    Y.clear();for(So t:y)Y.push_back(exp(t));
    heSo=binhPhuongToiThieu(x,Y,1); A=heSo[0]; B=heSo[1];
    a=B; b=A;
    cout << "Y=exp(y), X=x; Y=A+B*X; a=B, b=A\n";
    cout << "y=ln(a*x+b): a=" << a << "; b=" << b << '\n';
    inSaiSo("Ham logarit",x,y,[&](So t){return log(a*t+b);});
    cout << "Bac hai co S nho nhat trong bon dang da khop.\n";
}
void bai6() {
    cout << "\nBAI 6 - TAM THUC BAC HAI\n";
    Day x{.78L,1.56L,2.34L,3.12L,3.81L},y{2.50L,1.20L,1.12L,4.28L,7.50L};
    Day p=binhPhuongToiThieu(x,y,2);
    // Bai 6 dat y=a*x^2+b*x+c: a=p[2], b=p[1], c=p[0].
    cout << "y=a*x^2+b*x+c: a=" << p[2] << "; b=" << p[1] << "; c=" << p[0] << '\n';
    inSaiSo("Bac 2",x,y,[&](So t){return giaTri(p,t);});
}
void bai7() {
    cout << "\nBAI 7 - SAI SO NOI SUY SIN\n";
    So pi=acos(-1.0L);Day x,y;
    for(int d=0;d<=25;d+=5){So t=d*pi/180;x.push_back(t);y.push_back(sin(t));}
    cout << "x*=12 do 30 phut=12.5 do; tinh bang radian\n";
    So z=12.5L*pi/180,v=lagrangeTai(x,y,z),tich=1;
    for(So t:x)tich*=z-t;
    So chan=sin(25*pi/180)*abs(tich)/720;
    So e=abs(v-sin(z));
    cout << setprecision(15) << "p_5(x*)=" << v << "; sin=" << sin(z) << '\n';
    cout << scientific << "|R_5(x*)| quan sat=" << e << "; chan ly thuyet=" << chan << '\n'
              << fixed << setprecision(10);
    if(e>chan+1e-15L)throw runtime_error("Sai so vuot chan");
}
void bai8() {
    cout << "\nBAI 8 - SO SANH BA DANG\n";
    Day x{2,4,6,8,10,12},y{7.32L,8.24L,9.20L,10.19L,11.01L,12.05L};
    for(int d:{1,2}) {
        Day p=binhPhuongToiThieu(x,y,d);inMoHinh(p);
        inSaiSo(d==1?"Bac 1":"Bac 2",x,y,[&](So t){return giaTri(p,t);});
    }
    // Dung lg (logarit co so 10) dung nhu slide chuong 3, trang 13.
    Day X,Y;
    for(So t:x) X.push_back(log10(t));
    for(So t:y) Y.push_back(log10(t));
    Day heSo=binhPhuongToiThieu(X,Y,1);
    So A=heSo[0], B=heSo[1];
    So a=pow(10.0L,A), b=B;
    cout << "Y=lg(y), X=lg(x); Y=A+B*X\n";
    cout << "A=" << A << "; B=" << B << "; a=10^A, b=B\n";
    cout << "y=a*x^b: a=" << a << "; b=" << b << '\n';
    inSaiSo("Luy thua",x,y,[&](So t){return a*pow(t,b);});
    cout << "Bac 2 co S nho nhat; bac 1 gan tuong duong va don gian hon.\n";
}
int main() {
    cout << fixed << setprecision(10);
    cout << "BAI 9 SEMINAR 2 - KIEM TRA CAC BAI 1 DEN 8\n";
    cout << "Ky hieu: p_n, l_i (Lagrange); L_n (Newton); epsilon_i=y_i-f(x_i).\n";
    cout << "S=tong epsilon_i^2; lg=log10; ln=log tu nhien.\n";
    cout << "Cac dang phi tuyen duoc tuyen tinh hoa; S danh gia tren y goc.\n";
    try {bai1();bai2();bai3();bai4();bai5();bai6();bai7();bai8();}
    catch(const exception& e){cerr << "Loi: " << e.what() << '\n';return 1;}
    cout << "\nHoan thanh 8 bai.\n";
    return 0;
}
