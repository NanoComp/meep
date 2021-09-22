namespace meep {

// like custom_src_time, but using Python function object, with proper reference counting
class custom_py_src_time : public src_time {
public:
  custom_py_src_time(PyObject *fun, double st = -infinity, double et = infinity,
                     std::complex<double> f = 0)
      : func(fun), freq(f), start_time(float(st)), end_time(float(et)) {
    SWIG_PYTHON_THREAD_BEGIN_BLOCK;
    Py_INCREF(func);
    SWIG_PYTHON_THREAD_END_BLOCK;
  }
  virtual ~custom_py_src_time() {
    SWIG_PYTHON_THREAD_BEGIN_BLOCK;
    Py_DECREF(func);
    SWIG_PYTHON_THREAD_END_BLOCK;
  }

  virtual std::complex<double> current(double time, double dt) const {
    if (is_integrated)
      return src_time::current(time, dt);
    else
      return dipole(time);
  }
  virtual std::complex<double> dipole(double time) const {
    float rtime = float(time);
    if (rtime >= start_time && rtime <= end_time) {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      PyObject *py_t = PyFloat_FromDouble(time);
      PyObject *pyres = PyObject_CallFunctionObjArgs(func, py_t, NULL);
      double real = PyComplex_RealAsDouble(pyres);
      double imag = PyComplex_ImagAsDouble(pyres);
      std::complex<double> ret(real, imag);
      Py_DECREF(py_t);
      Py_DECREF(pyres);
      SWIG_PYTHON_THREAD_END_BLOCK;
      return ret;
    }
    else
      return 0.0;
  }
  virtual double last_time() const { return end_time; };
  virtual src_time *clone() const {
    SWIG_PYTHON_THREAD_BEGIN_BLOCK;
    Py_INCREF(func); // default copy constructor doesn't incref
    SWIG_PYTHON_THREAD_END_BLOCK;
    return new custom_py_src_time(*this);
  }
  virtual bool is_equal(const src_time &t) const {
    const custom_py_src_time *tp = dynamic_cast<const custom_py_src_time *>(&t);
    if (tp)
      return (tp->start_time == start_time && tp->end_time == end_time && tp->func == func &&
              tp->freq == freq);
    else
      return 0;
  }
  virtual std::complex<double> frequency() const { return freq; }
  virtual void set_frequency(std::complex<double> f) { freq = f; }

private:
  PyObject *func;
  std::complex<double> freq;
  double start_time, end_time;
};

} // namespace meep
