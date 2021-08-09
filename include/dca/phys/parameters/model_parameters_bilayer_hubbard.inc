// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for square lattice bilayer Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>> {
public:
  ModelParameters() : e1_(0.), e2_(0.), t1_(0.), t2_(0.), t1_prime_(0.), t2_prime_(0.), t_perp_(0.), U1_(0.), U2_(0.), V_(0.), V_prime_(0.) {}

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

  double get_t_perp() const {
    return t_perp_;
  }
  void set_t_perp(const double t_perp) {
    t_perp_ = t_perp;
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
  double t1_;
  double t2_;
  double t1_prime_;
  double t2_prime_;
  double t_perp_;
  double U1_;
  double U2_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(e1_);
  buffer_size += concurrency.get_buffer_size(e2_);
  buffer_size += concurrency.get_buffer_size(t1_);
  buffer_size += concurrency.get_buffer_size(t2_);
  buffer_size += concurrency.get_buffer_size(t1_prime_);
  buffer_size += concurrency.get_buffer_size(t2_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(U1_);
  buffer_size += concurrency.get_buffer_size(U2_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, e1_);
  concurrency.pack(buffer, buffer_size, position, e2_);
  concurrency.pack(buffer, buffer_size, position, t1_);
  concurrency.pack(buffer, buffer_size, position, t2_);
  concurrency.pack(buffer, buffer_size, position, t1_prime_);
  concurrency.pack(buffer, buffer_size, position, t2_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, U1_);
  concurrency.pack(buffer, buffer_size, position, U2_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, e1_);
  concurrency.unpack(buffer, buffer_size, position, e2_);
  concurrency.unpack(buffer, buffer_size, position, t1_);
  concurrency.unpack(buffer, buffer_size, position, t2_);
  concurrency.unpack(buffer, buffer_size, position, t1_prime_);
  concurrency.unpack(buffer, buffer_size, position, t2_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, U1_);
  concurrency.unpack(buffer, buffer_size, position, U2_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::bilayer_lattice<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("bilayer-Hubbard-model");

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
      reader_or_writer.execute("t-perp", t_perp_);
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
