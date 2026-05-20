// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang (May. 19, 2026)
// Trilayer two-orbital model for La-4310
//
// Trilayer lattice with totally 6 orbitals.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::La4Ni3O10<PointGroup>>> {
public:
  ModelParameters() : ex1_(0.), ez1_(0.), ex2_(0.), ez2_(0.), t1x1x1_(0.), t2x1x1_(0.), t3x1x1_(0.), t1x2x2_(0.), t2x2x2_(0.), t3x2x2_(0.), t1z1z1_(0.), t2z1z1_(0.), t3z1z1_(0.), t1z2z2_(0.), t2z2z2_(0.), t3z2z2_(0.), t1x1z1_(0.), t2x1z1_(0.), t1x2z2_(0.), t2x2z2_(0.), s0x1x2_(0.), s1x1x2_(0.), s2x1x2_(0.), s3x1x2_(0.), s0x2x2_(0.), s1x2x2_(0.), s2x2x2_(0.), s3x2x2_(0.), s0z1z2_(0.), s1z1z2_(0.), s2z1z2_(0.), s3z1z2_(0.), s0z2z2_(0.), s1z2z2_(0.), s2z2z2_(0.), s3z2z2_(0.), s1x1z2_(0.), s2x1z2_(0.), s1x2z2_(0.), s2x2z2_(0.), Ux_(0.), Uz_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_ex1() const {
    return ex1_;
  }
  void set_ex1(const double ex1) {
    ex1_ = ex1;
  }

  double get_ez1() const {
    return ez1_;
  }
  void set_ez1(const double ez1) {
    ez1_ = ez1;
  }

  double get_ex2() const {
    return ex2_;
  }
  void set_ex2(const double ex2) {
    ex2_ = ex2;
  }

  double get_ez2() const {
    return ez2_;
  }
  void set_ez2(const double ez2) {
    ez2_ = ez2;
  }

  double get_t1x1x1() const {
    return t1x1x1_;
  }
  void set_t1x1x1(const double t1x1x1) {
    t1x1x1_ = t1x1x1;
  }

  double get_t2x1x1() const {
    return t2x1x1_;
  }
  void set_t2x1x1(const double t2x1x1) {
    t2x1x1_ = t2x1x1;
  }

  double get_t3x1x1() const {
    return t3x1x1_;
  }
  void set_t3x1x1(const double t3x1x1) {
    t3x1x1_ = t3x1x1;
  }

  double get_t1x2x2() const {
    return t1x2x2_;
  }
  void set_t1x2x2(const double t1x2x2) {
    t1x2x2_ = t1x2x2;
  }

  double get_t2x2x2() const {
    return t2x2x2_;
  }
  void set_t2x2x2(const double t2x2x2) {
    t2x2x2_ = t2x2x2;
  }

  double get_t3x2x2() const {
    return t3x2x2_;
  }
  void set_t3x2x2(const double t3x2x2) {
    t3x2x2_ = t3x2x2;
  }

  double get_t1z1z1() const {
    return t1z1z1_;
  }
  void set_t1z1z1(const double t1z1z1) {
    t1z1z1_ = t1z1z1;
  }

  double get_t2z1z1() const {
    return t2z1z1_;
  }
  void set_t2z1z1(const double t2z1z1) {
    t2z1z1_ = t2z1z1;
  }

  double get_t3z1z1() const {
    return t3z1z1_;
  }
  void set_t3z1z1(const double t3z1z1) {
    t3z1z1_ = t3z1z1;
  }
  
   double get_t1z2z2() const {
    return t1z2z2_;
  }
  void set_t1z2z2(const double t1z2z2) {
    t1z2z2_ = t1z2z2;
  }

  double get_t2z2z2() const {
    return t2z2z2_;
  }
  void set_t2z2z2(const double t2z2z2) {
    t2z2z2_ = t2z2z2;
  }

  double get_t3z2z2() const {
    return t3z2z2_;
  }
  void set_t3z2z2(const double t3z2z2) {
    t3z2z2_ = t3z2z2;
  }
  
   double get_t1x1z1() const {
    return t1x1z1_;
  }
  void set_t1x1z1(const double t1x1z1) {
    t1x1z1_ = t1x1z1;
  }

  double get_t2x1z1() const {
    return t2x1z1_;
  }
  void set_t2x1z1(const double t2x1z1) {
    t2x1z1_ = t2x1z1;
  }

  double get_t1x2z2() const {
    return t1x2z2_;
  }
  void set_t1x2z2(const double t1x2z2) {
    t1x2z2_ = t1x2z2;
  }

  double get_t2x2z2() const {
    return t2x2z2_;
  }
  void set_t2x2z2(const double t2x2z2) {
    t2x2z2_ = t2x2z2;
  }

  double get_s0x1x2() const {
    return s0x1x2_;
  }
  void set_s0x1x2(const double s0x1x2) {
    s0x1x2_ = s0x1x2;
  }

  double get_s1x1x2() const {
    return s1x1x2_;
  }
  void set_s1x1x2(const double s1x1x2) {
    s1x1x2_ = s1x1x2;
  }

  double get_s2x1x2() const {
    return s2x1x2_;
  }
  void set_s2x1x2(const double s2x1x2) {
    s2x1x2_ = s2x1x2;
  }

  double get_s3x1x2() const {
    return s3x1x2_;
  }
  void set_s3x1x2(const double s3x1x2) {
    s3x1x2_ = s3x1x2;
  }

  double get_s0x2x2() const {
    return s0x2x2_;
  }
  void set_s0x2x2(const double s0x2x2) {
    s0x2x2_ = s0x2x2;
  }

  double get_s1x2x2() const {
    return s1x2x2_;
  }
  void set_s1x2x2(const double s1x2x2) {
    s1x2x2_ = s1x2x2;
  }

  double get_s2x2x2() const {
    return s2x2x2_;
  }
  void set_s2x2x2(const double s2x2x2) {
    s2x2x2_ = s2x2x2;
  }

  double get_s3x2x2() const {
    return s3x2x2_;
  }
  void set_s3x2x2(const double s3x2x2) {
    s3x2x2_ = s3x2x2;
  }

  double get_s0z1z2() const {
    return s0z1z2_;
  }
  void set_s0z1z2(const double s0z1z2) {
    s0z1z2_ = s0z1z2;
  }

  double get_s1z1z2() const {
    return s1z1z2_;
  }
  void set_s1z1z2(const double s1z1z2) {
    s1z1z2_ = s1z1z2;
  }

  double get_s2z1z2() const {
    return s2z1z2_;
  }
  void set_s2z1z2(const double s2z1z2) {
    s2z1z2_ = s2z1z2;
  }
  
  double get_s3z1z2() const {
    return s3z1z2_;
  }
  void set_s3z1z2(const double s3z1z2) {
    s3z1z2_ = s3z1z2;
  }

  double get_s0z2z2() const {
    return s0z2z2_;
  }
  void set_s0z2z2(const double s0z2z2) {
    s0z2z2_ = s0z2z2;
  }

  double get_s1z2z2() const {
    return s1z2z2_;
  }
  void set_s1z2z2(const double s1z2z2) {
    s1z2z2_ = s1z2z2;
  }

  double get_s2z2z2() const {
    return s2z2z2_;
  }
  void set_s2z2z2(const double s2z2z2) {
    s2z2z2_ = s2z2z2;
  }

  double get_s3z2z2() const {
    return s3z2z2_;
  }
  void set_s3z2z2(const double s3z2z2) {
    s3z2z2_ = s3z2z2;
  }

  double get_s1x1z2() const {
    return s1x1z2_;
  }
  void set_s1x1z2(const double s1x1z2) {
    s1x1z2_ = s1x1z2;
  }

  double get_s2x1z2() const {
    return s2x1z2_;
  }
  void set_s2x1z2(const double s2x1z2) {
    s2x1z2_ = s2x1z2;
  }

  double get_s1x2z2() const {
    return s1x2z2_;
  }
  void set_s1x2z2(const double s1x2z2) {
    s1x2z2_ = s1x2z2;
  }

  double get_s2x2z2() const {
    return s2x2z2_;
  }
  void set_s2x2z2(const double s2x2z2) {
    s2x2z2_ = s2x2z2;
  }

  double get_Ux() const {
    return Ux_;
  }
  void set_Ux(const double Ux) {
    Ux_ = Ux;
  }

  double get_Uz() const {
    return Uz_;
  }
  void set_Uz(const double Uz) {
    Uz_ = Uz;
  }

  double get_V() const {
    return V_;
  }
  void set_V(const double V) {
    V_ = V;
  }

  double get_V_prime() const {
    return V_prime_;
  }
  void set_V_prime(const double V_prime) {
    V_prime_ = V_prime;
  }

private:
  double ex1_;
  double ez1_;
  double ex2_;
  double ez2_;

  double t1x1x1_;
  double t2x1x1_;
  double t3x1x1_;
  double t1x2x2_;
  double t2x2x2_;
  double t3x2x2_;

  double t1z1z1_;
  double t2z1z1_;
  double t3z1z1_;
  double t1z2z2_;
  double t2z2z2_;
  double t3z2z2_;

  double t1x1z1_;
  double t2x1z1_;
  double t1x2z2_;
  double t2x2z2_;

  double s0x1x2_;
  double s1x1x2_;
  double s2x1x2_;
  double s3x1x2_;
  double s0x2x2_;
  double s1x2x2_;
  double s2x2x2_;
  double s3x2x2_;

  double s0z1z2_;
  double s1z1z2_;
  double s2z1z2_;
  double s3z1z2_;
  double s0z2z2_;
  double s1z2z2_;
  double s2z2z2_;
  double s3z2z2_;

  double s1x1z2_;
  double s2x1z2_;
  double s1x2z2_;
  double s2x2z2_;
  double Ux_;
  double Uz_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::La4Ni3O10<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ex1_);
  buffer_size += concurrency.get_buffer_size(ez1_);
  buffer_size += concurrency.get_buffer_size(ex2_);
  buffer_size += concurrency.get_buffer_size(ez2_);

  buffer_size += concurrency.get_buffer_size(t1x1x1_);
  buffer_size += concurrency.get_buffer_size(t2x1x1_);
  buffer_size += concurrency.get_buffer_size(t3x1x1_);
  buffer_size += concurrency.get_buffer_size(t1x2x2_);
  buffer_size += concurrency.get_buffer_size(t2x2x2_);
  buffer_size += concurrency.get_buffer_size(t3x2x2_);

  buffer_size += concurrency.get_buffer_size(t1z1z1_);
  buffer_size += concurrency.get_buffer_size(t2z1z1_);
  buffer_size += concurrency.get_buffer_size(t3z1z1_);
  buffer_size += concurrency.get_buffer_size(t1z2z2_);
  buffer_size += concurrency.get_buffer_size(t2z2z2_);
  buffer_size += concurrency.get_buffer_size(t3z2z2_);

  buffer_size += concurrency.get_buffer_size(t1x1z1_);
  buffer_size += concurrency.get_buffer_size(t2x1z1_);
  buffer_size += concurrency.get_buffer_size(t1x2z2_);
  buffer_size += concurrency.get_buffer_size(t2x2z2_);

  buffer_size += concurrency.get_buffer_size(s0x1x2_);
  buffer_size += concurrency.get_buffer_size(s1x1x2_);
  buffer_size += concurrency.get_buffer_size(s2x1x2_);
  buffer_size += concurrency.get_buffer_size(s3x1x2_);
  buffer_size += concurrency.get_buffer_size(s0x2x2_);
  buffer_size += concurrency.get_buffer_size(s1x2x2_);
  buffer_size += concurrency.get_buffer_size(s2x2x2_);
  buffer_size += concurrency.get_buffer_size(s3x2x2_);

  buffer_size += concurrency.get_buffer_size(s0z1z2_);
  buffer_size += concurrency.get_buffer_size(s1z1z2_);
  buffer_size += concurrency.get_buffer_size(s2z1z2_);
  buffer_size += concurrency.get_buffer_size(s3z1z2_);
  buffer_size += concurrency.get_buffer_size(s0z2z2_);
  buffer_size += concurrency.get_buffer_size(s1z2z2_);
  buffer_size += concurrency.get_buffer_size(s2z2z2_);
  buffer_size += concurrency.get_buffer_size(s3z2z2_);
 
  buffer_size += concurrency.get_buffer_size(s1x1z2_);
  buffer_size += concurrency.get_buffer_size(s2x1z2_);
  buffer_size += concurrency.get_buffer_size(s1x2z2_);
  buffer_size += concurrency.get_buffer_size(s2x2z2_);
  buffer_size += concurrency.get_buffer_size(Ux_);
  buffer_size += concurrency.get_buffer_size(Uz_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::La4Ni3O10<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, ex1_);
  concurrency.pack(buffer, buffer_size, position, ez1_);
  concurrency.pack(buffer, buffer_size, position, ex2_);
  concurrency.pack(buffer, buffer_size, position, ez2_);

  concurrency.pack(buffer, buffer_size, position, t1x1x1_);
  concurrency.pack(buffer, buffer_size, position, t2x1x1_);
  concurrency.pack(buffer, buffer_size, position, t3x1x1_);
  concurrency.pack(buffer, buffer_size, position, t1x2x2_);
  concurrency.pack(buffer, buffer_size, position, t2x2x2_);
  concurrency.pack(buffer, buffer_size, position, t3x2x2_);
  
  concurrency.pack(buffer, buffer_size, position, t1z1z1_);
  concurrency.pack(buffer, buffer_size, position, t2z1z1_);
  concurrency.pack(buffer, buffer_size, position, t3z1z1_);
  concurrency.pack(buffer, buffer_size, position, t1z2z2_);
  concurrency.pack(buffer, buffer_size, position, t2z2z2_);
  concurrency.pack(buffer, buffer_size, position, t3z2z2_);

  concurrency.pack(buffer, buffer_size, position, t1x1z1_);
  concurrency.pack(buffer, buffer_size, position, t2x1z1_);
  concurrency.pack(buffer, buffer_size, position, t1x2z2_);
  concurrency.pack(buffer, buffer_size, position, t2x2z2_);

  concurrency.pack(buffer, buffer_size, position, s0x1x2_);
  concurrency.pack(buffer, buffer_size, position, s1x1x2_);
  concurrency.pack(buffer, buffer_size, position, s2x1x2_);
  concurrency.pack(buffer, buffer_size, position, s3x1x2_);
  concurrency.pack(buffer, buffer_size, position, s0x2x2_);
  concurrency.pack(buffer, buffer_size, position, s1x2x2_);
  concurrency.pack(buffer, buffer_size, position, s2x2x2_);
  concurrency.pack(buffer, buffer_size, position, s3x2x2_);

  concurrency.pack(buffer, buffer_size, position, s0z1z2_);
  concurrency.pack(buffer, buffer_size, position, s1z1z2_);
  concurrency.pack(buffer, buffer_size, position, s2z1z2_);
  concurrency.pack(buffer, buffer_size, position, s3z1z2_);
  concurrency.pack(buffer, buffer_size, position, s0z2z2_);
  concurrency.pack(buffer, buffer_size, position, s1z2z2_);
  concurrency.pack(buffer, buffer_size, position, s2z2z2_);
  concurrency.pack(buffer, buffer_size, position, s3z2z2_);

  concurrency.pack(buffer, buffer_size, position, s1x1z2_);
  concurrency.pack(buffer, buffer_size, position, s2x1z2_);
  concurrency.pack(buffer, buffer_size, position, s1x2z2_);
  concurrency.pack(buffer, buffer_size, position, s2x2z2_);
  concurrency.pack(buffer, buffer_size, position, Ux_);
  concurrency.pack(buffer, buffer_size, position, Uz_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::La4Ni3O10<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, ex1_);
  concurrency.unpack(buffer, buffer_size, position, ez1_);
  concurrency.unpack(buffer, buffer_size, position, ex2_);
  concurrency.unpack(buffer, buffer_size, position, ez2_);

  concurrency.unpack(buffer, buffer_size, position, t1x1x1_);
  concurrency.unpack(buffer, buffer_size, position, t2x1x1_);
  concurrency.unpack(buffer, buffer_size, position, t3x1x1_);
  concurrency.unpack(buffer, buffer_size, position, t1x2x2_);
  concurrency.unpack(buffer, buffer_size, position, t2x2x2_);
  concurrency.unpack(buffer, buffer_size, position, t3x2x2_);

  concurrency.unpack(buffer, buffer_size, position, t1z1z1_);
  concurrency.unpack(buffer, buffer_size, position, t2z1z1_);
  concurrency.unpack(buffer, buffer_size, position, t3z1z1_);
  concurrency.unpack(buffer, buffer_size, position, t1z2z2_);
  concurrency.unpack(buffer, buffer_size, position, t2z2z2_);
  concurrency.unpack(buffer, buffer_size, position, t3z2z2_);

  concurrency.unpack(buffer, buffer_size, position, t1x1z1_);
  concurrency.unpack(buffer, buffer_size, position, t2x1z1_);
  concurrency.unpack(buffer, buffer_size, position, t1x2z2_);
  concurrency.unpack(buffer, buffer_size, position, t2x2z2_);

  concurrency.unpack(buffer, buffer_size, position, s0x1x2_);
  concurrency.unpack(buffer, buffer_size, position, s1x1x2_);
  concurrency.unpack(buffer, buffer_size, position, s2x1x2_);
  concurrency.unpack(buffer, buffer_size, position, s3x1x2_);
  concurrency.unpack(buffer, buffer_size, position, s0x2x2_);
  concurrency.unpack(buffer, buffer_size, position, s1x2x2_);
  concurrency.unpack(buffer, buffer_size, position, s2x2x2_);
  concurrency.unpack(buffer, buffer_size, position, s3x2x2_);

  concurrency.unpack(buffer, buffer_size, position, s0z1z2_);
  concurrency.unpack(buffer, buffer_size, position, s1z1z2_);
  concurrency.unpack(buffer, buffer_size, position, s2z1z2_);
  concurrency.unpack(buffer, buffer_size, position, s3z1z2_);
  concurrency.unpack(buffer, buffer_size, position, s0z2z2_);
  concurrency.unpack(buffer, buffer_size, position, s1z2z2_);
  concurrency.unpack(buffer, buffer_size, position, s2z2z2_);
  concurrency.unpack(buffer, buffer_size, position, s3z2z2_);

  concurrency.unpack(buffer, buffer_size, position, s1x1z2_);
  concurrency.unpack(buffer, buffer_size, position, s2x1z2_);
  concurrency.unpack(buffer, buffer_size, position, s1x2z2_);
  concurrency.unpack(buffer, buffer_size, position, s2x2z2_);
  concurrency.unpack(buffer, buffer_size, position, Ux_);
  concurrency.unpack(buffer, buffer_size, position, Uz_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::La4Ni3O10<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("La4Ni3O10-model");

    try {
      reader_or_writer.execute("ex1", ex1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ez1", ez1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ex2", ex2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ez2", ez2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1x1x1", t1x1x1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2x1x1", t2x1x1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3x1x1", t3x1x1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1x2x2", t1x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2x2x2", t2x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3x2x2", t3x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1z1z1", t1z1z1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2z1z1", t2z1z1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3z1z1", t3z1z1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1z2z2", t1z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2z2z2", t2z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3z2z2", t3z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1x1z1", t1x1z1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2x1z1", t2x1z1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1x2z2", t1x2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2x2z2", t2x2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s0x1x2", s0x1x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1x1x2", s1x1x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2x1x2", s2x1x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s3x1x2", s3x1x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s0x2x2", s0x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1x2x2", s1x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2x2x2", s2x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s3x2x2", s3x2x2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s0z1z2", s0z1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1z1z2", s1z1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2z1z2", s2z1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s3z1z2", s3z1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s0z2z2", s0z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1z2z2", s1z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2z2z2", s2z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s3z2z2", s3z2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1x1z2", s1x1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2x1z2", s2x1z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s1x2z2", s1x2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s2x2z2", s2x2z2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("Ux", Ux_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("Uz", Uz_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V", V_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V-prime", V_prime_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
