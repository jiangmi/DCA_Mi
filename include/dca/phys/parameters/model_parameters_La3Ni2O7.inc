// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang (Aug. 29, 2025)
// Bilayer two-orbital model for La-327
//
// Bilayer lattice with totally 4 orbitals.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::La3Ni2O7<PointGroup>>> {
public:
  ModelParameters() : ex_(0.), ez_(0.), t11x_(0.), t11xy_(0.), t11xx_(0.), t22x_(0.), t22xy_(0.), t22xx_(0.), t12x_(0.), t12xx_(0.), s110_(0.), s11x_(0.), s11xy_(0.), s11xx_(0.), s220_(0.), s22x_(0.), s22xy_(0.), s22xx_(0.), s12x_(0.), s12xx_(0.), Ux_(0.), Uz_(0.), V_(0.), V_prime_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_ex() const {
    return ex_;
  }
  void set_ex(const double ex) {
    ex_ = ex;
  }

  double get_ez() const {
    return ez_;
  }
  void set_ez(const double ez) {
    ez_ = ez;
  }

  double get_t11x() const {
    return t11x_;
  }
  void set_t11x(const double t11x) {
    t11x_ = t11x;
  }

  double get_t11xy() const {
    return t11xy_;
  }
  void set_t11xy(const double t11xy) {
    t11xy_ = t11xy;
  }
  
  double get_t11xx() const {
    return t11xx_;
  }
  void set_t11xx(const double t11xx) {
    t11xx_ = t11xx;
  }
  
  double get_t22x() const {
    return t22x_;
  }
  void set_t22x(const double t22x) {
    t22x_ = t22x;
  }

  double get_t22xy() const {
    return t22xy_;
  }
  void set_t22xy(const double t22xy) {
    t22xy_ = t22xy;
  }
  
  double get_t22xx() const {
    return t22xx_;
  }
  void set_t22xx(const double t22xx) {
    t22xx_ = t22xx;
  }
  
  double get_t12x() const {
    return t12x_;
  }
  void set_t12x(const double t12x) {
    t12x_ = t12x;
  }
  
  double get_t12xx() const {
    return t12xx_;
  }
  void set_t12xx(const double t12xx) {
    t12xx_ = t12xx;
  }
  
  double get_s110() const {
    return s110_;
  }
  void set_s110(const double s110) {
    s110_ = s110;
  }
  
  double get_s11x() const {
    return s11x_;
  }
  void set_s11x(const double s11x) {
    s11x_ = s11x;
  }

  double get_s11xy() const {
    return s11xy_;
  }
  void set_s11xy(const double s11xy) {
    s11xy_ = s11xy;
  }
  
  double get_s11xx() const {
    return s11xx_;
  }
  void set_s11xx(const double s11xx) {
    s11xx_ = s11xx;
  }
  
  double get_s220() const {
    return s220_;
  }
  void set_s220(const double s220) {
    s220_ = s220;
  }
  
  double get_s22x() const {
    return s22x_;
  }
  void set_s22x(const double s22x) {
    s22x_ = s22x;
  }

  double get_s22xy() const {
    return s22xy_;
  }
  void set_s22xy(const double s22xy) {
    s22xy_ = s22xy;
  }
  
  double get_s22xx() const {
    return s22xx_;
  }
  void set_s22xx(const double s22xx) {
    s22xx_ = s22xx;
  }
  
  double get_s12x() const {
    return s12x_;
  }
  void set_s12x(const double s12x) {
    s12x_ = s12x;
  }
  
  double get_s12xx() const {
    return s12xx_;
  }
  void set_s12xx(const double s12xx) {
    s12xx_ = s12xx;
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
  double ex_;
  double ez_;
  double t11x_;
  double t11xy_;
  double t11xx_;
  double t22x_;
  double t22xy_;
  double t22xx_;
  double t12x_;
  double t12xx_;
  double s110_; 
  double s11x_;
  double s11xy_;
  double s11xx_;
  double s220_; 
  double s22x_;
  double s22xy_;
  double s22xx_;
  double s12x_;
  double s12xx_;
  double Ux_;
  double Uz_;
  double V_;
  double V_prime_;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::La3Ni2O7<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(ex_);
  buffer_size += concurrency.get_buffer_size(ez_);
  buffer_size += concurrency.get_buffer_size(t11x_);
  buffer_size += concurrency.get_buffer_size(t11xy_);
  buffer_size += concurrency.get_buffer_size(t11xx_);
  buffer_size += concurrency.get_buffer_size(t22x_);
  buffer_size += concurrency.get_buffer_size(t22xy_);
  buffer_size += concurrency.get_buffer_size(t22xx_);
  buffer_size += concurrency.get_buffer_size(t12x_);
  buffer_size += concurrency.get_buffer_size(t12xx_);
  buffer_size += concurrency.get_buffer_size(s110_);
  buffer_size += concurrency.get_buffer_size(s11x_);
  buffer_size += concurrency.get_buffer_size(s11xy_);
  buffer_size += concurrency.get_buffer_size(s11xx_);
  buffer_size += concurrency.get_buffer_size(s220_);
  buffer_size += concurrency.get_buffer_size(s22x_);
  buffer_size += concurrency.get_buffer_size(s22xy_);
  buffer_size += concurrency.get_buffer_size(s22xx_);
  buffer_size += concurrency.get_buffer_size(s12x_);
  buffer_size += concurrency.get_buffer_size(s12xx_);
  buffer_size += concurrency.get_buffer_size(Ux_);
  buffer_size += concurrency.get_buffer_size(Uz_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::La3Ni2O7<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, ex_);
  concurrency.pack(buffer, buffer_size, position, ez_);
  concurrency.pack(buffer, buffer_size, position, t11x_);
  concurrency.pack(buffer, buffer_size, position, t11xy_);
  concurrency.pack(buffer, buffer_size, position, t11xx_);
  concurrency.pack(buffer, buffer_size, position, t22x_);
  concurrency.pack(buffer, buffer_size, position, t22xy_);
  concurrency.pack(buffer, buffer_size, position, t22xx_);
  concurrency.pack(buffer, buffer_size, position, t12x_);
  concurrency.pack(buffer, buffer_size, position, t12xx_);
  concurrency.pack(buffer, buffer_size, position, s110_);
  concurrency.pack(buffer, buffer_size, position, s11x_);
  concurrency.pack(buffer, buffer_size, position, s11xy_);
  concurrency.pack(buffer, buffer_size, position, s11xx_);
  concurrency.pack(buffer, buffer_size, position, s220_);
  concurrency.pack(buffer, buffer_size, position, s22x_);
  concurrency.pack(buffer, buffer_size, position, s22xy_);
  concurrency.pack(buffer, buffer_size, position, s22xx_);
  concurrency.pack(buffer, buffer_size, position, s12x_);
  concurrency.pack(buffer, buffer_size, position, s12xx_);
  concurrency.pack(buffer, buffer_size, position, Ux_);
  concurrency.pack(buffer, buffer_size, position, Uz_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::La3Ni2O7<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, ex_);
  concurrency.unpack(buffer, buffer_size, position, ez_);
  concurrency.unpack(buffer, buffer_size, position, t11x_);
  concurrency.unpack(buffer, buffer_size, position, t11xy_);
  concurrency.unpack(buffer, buffer_size, position, t11xx_);
  concurrency.unpack(buffer, buffer_size, position, t22x_);
  concurrency.unpack(buffer, buffer_size, position, t22xy_);
  concurrency.unpack(buffer, buffer_size, position, t22xx_);
  concurrency.unpack(buffer, buffer_size, position, t12x_);
  concurrency.unpack(buffer, buffer_size, position, t12xx_);
  concurrency.unpack(buffer, buffer_size, position, s110_);
  concurrency.unpack(buffer, buffer_size, position, s11x_);
  concurrency.unpack(buffer, buffer_size, position, s11xy_);
  concurrency.unpack(buffer, buffer_size, position, s11xx_);
  concurrency.unpack(buffer, buffer_size, position, s220_);
  concurrency.unpack(buffer, buffer_size, position, s22x_);
  concurrency.unpack(buffer, buffer_size, position, s22xy_);
  concurrency.unpack(buffer, buffer_size, position, s22xx_);
  concurrency.unpack(buffer, buffer_size, position, s12x_);
  concurrency.unpack(buffer, buffer_size, position, s12xx_);
  concurrency.unpack(buffer, buffer_size, position, Ux_);
  concurrency.unpack(buffer, buffer_size, position, Uz_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::La3Ni2O7<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("La3Ni2O7-model");

    try {
      reader_or_writer.execute("ex", ex_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ez", ez_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t11x", t11x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t11xy", t11xy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t11xx", t11xx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t22x", t22x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t22xy", t22xy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t22xx", t22xx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t12x", t12x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t12xx", t12xx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s110", s110_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s11x", s11x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s11xy", s11xy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s11xx", s11xx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s220", s220_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s22x", s22x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s22xy", s22xy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s22xx", s22xx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s12x", s12x_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("s12xx", s12xx_);
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
