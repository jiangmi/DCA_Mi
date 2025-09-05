// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang
//
// Template specialization of model parameters for tetralayer lattice.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::tetralayer_lattice<PointGroup>>> {
public:
  ModelParameters() : e1_(0.), e2_(0.), e3_(0.), e4_(0.), t1_(0.), t2_(0.), t3_(0.), t4_(0.), t1_prime_(0.), t2_prime_(0.), t3_prime_(0.), t4_prime_(0.), t_perp12_(0.), t_perp_prime12_(0.), t_perp23_(0.), t_perp_prime23_(0.), t_perp34_(0.), t_perp_prime34_(0.), coefp_(0.), coefm_(0.), U1_(0.), U2_(0.), U3_(0.), U4_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_e1() const {
    return e1_;
  }
  void set_e1(const double e1) {
    e1_ = e1;
  }

  double get_e2() const {
    return e2_;
  }
  void set_e2(const double e2) {
    e2_ = e2;
  }

  double get_e3() const {
    return e3_;
  }
  void set_e3(const double e3) {
    e3_ = e3;
  }

  double get_e4() const { 
    return e4_;
  }
  void set_e4(const double e4) {
   e4_ = e4;
  }

  double get_t1() const {
    return t1_;
  }
  void set_t1(const double t1) {
    t1_ = t1;
  }

  double get_t2() const {
    return t2_;
  }
  void set_t2(const double t2) {
    t2_ = t2;
  }

  double get_t3() const {
    return t3_;
  }
  void set_t3(const double t3) {
    t3_ = t3;
  }
  
   double get_t4() const {
    return t4_;
  }
   void set_t4(const double t4) {
    t4_ = t4;
  }

  double get_t1_prime() const {
    return t1_prime_;
  }
  void set_t1_prime(const double t1_prime) {
    t1_prime_ = t1_prime;
  }

  double get_t2_prime() const {
    return t2_prime_;
  }
  void set_t2_prime(const double t2_prime) {
    t2_prime_ = t2_prime;
  }

  double get_t3_prime() const {
    return t3_prime_;
  }
  void set_t3_prime(const double t3_prime) {
    t3_prime_ = t3_prime;
  }

  double get_t4_prime() const {
    return t4_prime_;
  }
  void set_t4_prime(const double t4_prime) {
   t4_prime_ = t4_prime;
  }

  double get_t_perp12() const {
    return t_perp12_;
  }
  void set_t_perp12(const double t_perp12) {
    t_perp12_ = t_perp12;
  }

  double get_t_perp_prime12() const {
    return t_perp_prime12_;
  }
  void set_t_perp_prime12(const double t_perp_prime12) {
    t_perp_prime12_ = t_perp_prime12;
  }
  
  double get_t_perp23() const {
    return t_perp23_;
  }
  void set_t_perp23(const double t_perp23) {
    t_perp23_ = t_perp23;
  }

  double get_t_perp_prime23() const {
    return t_perp_prime23_;
  }
  void set_t_perp_prime23(const double t_perp_prime23) {
    t_perp_prime23_ = t_perp_prime23;
  }

  double get_t_perp34() const {
    return t_perp34_;
  }
  void set_t_perp34(const double t_perp34) {
    t_perp34_ = t_perp34;
  }

  double get_t_perp_prime34() const {
    return t_perp_prime34_;
  }
  void set_t_perp_prime34(const double t_perp_prime34) {
    t_perp_prime34_ = t_perp_prime34;
 }
  
  double get_coefp() const {
    return coefp_;
  }
  void set_coefp(const double coefp) {
    coefp_ = coefp;
  }

  double get_coefm() const {
    return coefm_;
  }
  void set_coefm(const double coefm) {
    coefm_ = coefm;
  }

  double get_U1() const {
    return U1_;
  }
  void set_U1(const double U1) {
    U1_ = U1;
  }

  double get_U2() const {
    return U2_;
  }
  void set_U2(const double U2) {
    U2_ = U2;
  }

  double get_U3() const {
    return U3_;
  }
  void set_U3(const double U3) {
    U3_ = U3;
  }

  double get_U4() const {
    return U4_;
  }
  void set_U4(const double U4) {
    U4_ = U4;
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
  double e1_;
  double e2_;
  double e3_;
  double e4_;
  double t1_;
  double t2_;
  double t3_;
  double t4_;
  double t1_prime_;
  double t2_prime_;
  double t3_prime_;
  double t4_prime_;
  double t_perp12_;
  double t_perp_prime12_;
  double t_perp23_;
  double t_perp_prime23_;
  double t_perp34_;
  double t_perp_prime34_;
  double coefp_;
  double coefm_;
  double U1_;
  double U2_;
  double U3_;
  double U4_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::tetralayer_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(e1_);
  buffer_size += concurrency.get_buffer_size(e2_);
  buffer_size += concurrency.get_buffer_size(e3_);
  buffer_size += concurrency.get_buffer_size(e4_);
  buffer_size += concurrency.get_buffer_size(t1_);
  buffer_size += concurrency.get_buffer_size(t2_);
  buffer_size += concurrency.get_buffer_size(t3_);
  buffer_size += concurrency.get_buffer_size(t4_);
  buffer_size += concurrency.get_buffer_size(t1_prime_);
  buffer_size += concurrency.get_buffer_size(t2_prime_);
  buffer_size += concurrency.get_buffer_size(t3_prime_);
  buffer_size += concurrency.get_buffer_size(t4_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp12_);
  buffer_size += concurrency.get_buffer_size(t_perp_prime12_);
  buffer_size += concurrency.get_buffer_size(t_perp23_);
  buffer_size += concurrency.get_buffer_size(t_perp_prime23_);
  buffer_size += concurrency.get_buffer_size(t_perp34_);
  buffer_size += concurrency.get_buffer_size(t_perp_prime34_);
  buffer_size += concurrency.get_buffer_size(coefp_);
  buffer_size += concurrency.get_buffer_size(coefm_);
  buffer_size += concurrency.get_buffer_size(U1_);
  buffer_size += concurrency.get_buffer_size(U2_);
  buffer_size += concurrency.get_buffer_size(U3_);
  buffer_size += concurrency.get_buffer_size(U4_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::tetralayer_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, e1_);
  concurrency.pack(buffer, buffer_size, position, e2_);
  concurrency.pack(buffer, buffer_size, position, e3_);
  concurrency.pack(buffer, buffer_size, position, e4_);
  concurrency.pack(buffer, buffer_size, position, t1_);
  concurrency.pack(buffer, buffer_size, position, t2_);
  concurrency.pack(buffer, buffer_size, position, t3_);
  concurrency.pack(buffer, buffer_size, position, t4_);
  concurrency.pack(buffer, buffer_size, position, t1_prime_);
  concurrency.pack(buffer, buffer_size, position, t2_prime_);
  concurrency.pack(buffer, buffer_size, position, t3_prime_);
  concurrency.pack(buffer, buffer_size, position, t4_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp12_);
  concurrency.pack(buffer, buffer_size, position, t_perp_prime12_);
  concurrency.pack(buffer, buffer_size, position, t_perp23_);
  concurrency.pack(buffer, buffer_size, position, t_perp_prime23_);
  concurrency.pack(buffer, buffer_size, position, t_perp34_);
  concurrency.pack(buffer, buffer_size, position, t_perp_prime34_);
  concurrency.pack(buffer, buffer_size, position, coefp_);
  concurrency.pack(buffer, buffer_size, position, coefm_);
  concurrency.pack(buffer, buffer_size, position, U1_);
  concurrency.pack(buffer, buffer_size, position, U2_);
  concurrency.pack(buffer, buffer_size, position, U3_);
  concurrency.pack(buffer, buffer_size, position, U4_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::tetralayer_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, e1_);
  concurrency.unpack(buffer, buffer_size, position, e2_);
  concurrency.unpack(buffer, buffer_size, position, e3_);
  concurrency.unpack(buffer, buffer_size, position, e4_);
  concurrency.unpack(buffer, buffer_size, position, t1_);
  concurrency.unpack(buffer, buffer_size, position, t2_);
  concurrency.unpack(buffer, buffer_size, position, t3_);
  concurrency.unpack(buffer, buffer_size, position, t4_);
  concurrency.unpack(buffer, buffer_size, position, t1_prime_);
  concurrency.unpack(buffer, buffer_size, position, t2_prime_);
  concurrency.unpack(buffer, buffer_size, position, t3_prime_);
  concurrency.unpack(buffer, buffer_size, position, t4_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp12_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_prime12_);
  concurrency.unpack(buffer, buffer_size, position, t_perp23_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_prime23_);
  concurrency.unpack(buffer, buffer_size, position, t_perp34_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_prime34_);
  concurrency.unpack(buffer, buffer_size, position, coefp_);
  concurrency.unpack(buffer, buffer_size, position, coefm_);
  concurrency.unpack(buffer, buffer_size, position, U1_);
  concurrency.unpack(buffer, buffer_size, position, U2_);
  concurrency.unpack(buffer, buffer_size, position, U3_);
  concurrency.unpack(buffer, buffer_size, position, U4_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::tetralayer_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("tetralayer-Hubbard-model");

    try {
      reader_or_writer.execute("e1", e1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("e2", e2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("e3", e3_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("e4", e4_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1", t1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2", t2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3", t3_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t4", t4_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t1-prime", t1_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t2-prime", t2_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t3-prime", t3_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t4-prime", t4_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp12", t_perp12_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp-prime12", t_perp_prime12_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp23", t_perp23_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp-prime23", t_perp_prime23_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp34", t_perp34_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-perp-prime34", t_perp_prime34_);
    }
    catch (const std::exception& r_e) {
    }

    try {
      reader_or_writer.execute("coefp", coefp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("coefm", coefm_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U1", U1_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U2", U2_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U3", U3_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U4", U4_);
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
