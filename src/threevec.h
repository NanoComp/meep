#define FOR3(i) for (int i=0;i<3;i++)

class threevec {
public:
  threevec() {};
  threevec(int i) {
    FOR3(j) val[j] = (i==j)?1.0:0.0;
  };
  double operator*(const threevec &a) const {
    double out;
    FOR3(i) out += val[i]*a.val[i];
    return out;
  }
  threevec operator^(const threevec &a) const {
    threevec out;
    FOR3(i) out.val[i] = val[(i+1)%3]*a.val[(i+2)%3]
      - val[(i+2)%3]*a.val[(i+1)%3];
    return out;
  }
  threevec operator/=(double s) {
    FOR3(i) val[i] = (1.0/s)*val[i];
    return *this;
  }
  threevec operator/(double s) const {
    threevec out;
    FOR3(i) out.val[i] = (1.0/s)*val[i];
    return out;
  }
  threevec operator*(double s) const {
    threevec out;
    FOR3(i) out.val[i] = s*val[i];
    return out;
  }
  threevec operator-(const threevec &a) const {
    threevec out;
    FOR3(i) out.val[i] = val[i] - a.val[i];
    return out;
  };
  threevec operator+(const threevec &a) const {
    threevec out;
    FOR3(i) out.val[i] = val[i] + a.val[i];
    return out;
  };
  double val[3];
};

double abs(const threevec &v) {
  return sqrt(abs(v*v));
}

class tensor {
public:
  tensor() {};
  tensor(double v) {
    FOR3(i) row[i] = threevec(i)*v;
  };
  tensor(const threevec &v) { // Projection tensor.
    FOR3(i) FOR3(j) row[i].val[j] = v.val[i]*v.val[j];
  };
  ~tensor() { return; };

  tensor transpose() const {
    tensor out;
    FOR3(i) FOR3(j) out.row[i].val[j] = row[j].val[i];
    return out;
  }
  tensor operator*(const tensor &t) {
    tensor out;
    tensor ttrans = t.transpose();
    FOR3(i) out.row[i] = (*this)*ttrans.row[i];
    return out.transpose();
  }
  threevec operator*(const threevec &v) {
    threevec out;
    FOR3(i) out.val[i] = row[i]*v;
    return out;
  }
  tensor operator/(double s) const {
    tensor out;
    FOR3(i) out.row[i] = row[i]*(1.0/s);
    return out;
  }
  tensor operator*(double s) const {
    tensor out;
    FOR3(i) out.row[i] = row[i]*s;
    return out;
  }
  tensor operator-(const tensor &a) const {
    tensor out;
    FOR3(i) out.row[i] = row[i] - a.row[i];
    return out;
  };
  tensor operator+(const tensor &a) const {
    tensor out;
    FOR3(i) out.row[i] = row[i] + a.row[i];
    return out;
  };
  tensor operator+=(const tensor &a) {
    FOR3(i) row[i] = row[i] + a.row[i];
    return *this;
  };
  threevec row[3];
};

inline tensor transpose(const tensor &t) {
  return t.transpose();
}

inline double trace(const tensor &t) {
  double out = 0.0;
  FOR3(i) out += t.row[i].val[i];
  return out;
}

inline tensor operator/(double s, const tensor &t) {
  tensor out;
  double s_o_vol = s/(t.row[0]*(t.row[1]^t.row[2]));
  if (fabs(s_o_vol) > 1e50) {
    FOR3(i) out.row[i] = t.row[(i+1)%3]^t.row[(i+2)%3];
    FOR3(i) FOR3(j)
      if (out.row[i].val[j] != 0.0) out.row[i].val[j] = 1e80;
  } else {
    FOR3(i) out.row[i] = t.row[(i+1)%3]^t.row[(i+2)%3]*s_o_vol;
  }
  FOR3(i) FOR3(j)
    if (fabs(out.row[i].val[j]) < 1e-50) out.row[i].val[j] = 0.0;
  return out;
}

