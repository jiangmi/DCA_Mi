// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for the single-band Hubbard model.

template <typename Lattice>
class ModelParameters<models::TightBindingModel<Lattice>> {
public:
  ModelParameters() : tx_(0.), ty_(0.), tx_nn_(0.), ty_nn_(0.), t_prime_(0.), U_(0.), V_(0.), V_prime_(0.), Vpp_(0.) {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_tx() const {
    return tx_;
  }
  void set_tx(const double tx) {
    tx_ = tx;
  }

  double get_ty() const {
    return ty_;
  }
  void set_ty(const double ty) {
    ty_ = ty;
  }

  double get_tx_nn() const {
    return tx_nn_;
  }
  void set_tx_nn(const double tx_nn) {
    tx_nn_ = tx_nn;
  }

  double get_ty_nn() const {
    return ty_nn_;
  }
  void set_ty_nn(const double ty_nn) {
    ty_nn_ = ty_nn;
  }

  double get_t_prime() const {
    return t_prime_;
  }
  void set_t_prime(const double t_prime) {
    t_prime_ = t_prime;
  }

  double get_U() const {
    return U_;
  }
  void set_U(const double U) {
    U_ = U;
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

  double get_Vpp() const {
    return Vpp_;
  }
  void set_Vpp(const double Vpp) {
    Vpp_ = Vpp;
  }

private:
  double tx_;
  double ty_;
  double tx_nn_;
  double ty_nn_;
  double t_prime_;
  double U_;
  double V_;
  double V_prime_;
  double Vpp_;
};

template <typename Lattice>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<Lattice>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(tx_);
  buffer_size += concurrency.get_buffer_size(ty_);
  buffer_size += concurrency.get_buffer_size(tx_nn_);
  buffer_size += concurrency.get_buffer_size(ty_nn_);
  buffer_size += concurrency.get_buffer_size(t_prime_);
  buffer_size += concurrency.get_buffer_size(U_);
  buffer_size += concurrency.get_buffer_size(V_);
  buffer_size += concurrency.get_buffer_size(V_prime_);
  buffer_size += concurrency.get_buffer_size(Vpp_);

  return buffer_size;
}

template <typename Lattice>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<Lattice>>::pack(const Concurrency& concurrency,
                                                               char* buffer, int buffer_size,
                                                               int& position) const {
  concurrency.pack(buffer, buffer_size, position, tx_);
  concurrency.pack(buffer, buffer_size, position, ty_);
  concurrency.pack(buffer, buffer_size, position, tx_nn_);
  concurrency.pack(buffer, buffer_size, position, ty_nn_);
  concurrency.pack(buffer, buffer_size, position, t_prime_);
  concurrency.pack(buffer, buffer_size, position, U_);
  concurrency.pack(buffer, buffer_size, position, V_);
  concurrency.pack(buffer, buffer_size, position, V_prime_);
  concurrency.pack(buffer, buffer_size, position, Vpp_);
}

template <typename Lattice>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<Lattice>>::unpack(const Concurrency& concurrency,
                                                                 char* buffer, int buffer_size,
                                                                 int& position) {
  concurrency.unpack(buffer, buffer_size, position, tx_);
  concurrency.unpack(buffer, buffer_size, position, ty_);
  concurrency.unpack(buffer, buffer_size, position, tx_nn_);
  concurrency.unpack(buffer, buffer_size, position, ty_nn_);
  concurrency.unpack(buffer, buffer_size, position, t_prime_);
  concurrency.unpack(buffer, buffer_size, position, U_);
  concurrency.unpack(buffer, buffer_size, position, V_);
  concurrency.unpack(buffer, buffer_size, position, V_prime_);
  concurrency.unpack(buffer, buffer_size, position, Vpp_);
}

template <typename Lattice>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<Lattice>>::readWrite(ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("single-band-Hubbard-model");

    try {
      reader_or_writer.execute("tx", tx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ty", ty_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("tx_nn", tx_nn_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ty_nn", ty_nn_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t-prime", t_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U", U_);
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
    try {
      reader_or_writer.execute("Vpp", Vpp_);
    }
    catch (const std::exception& r_e) {
    }

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
