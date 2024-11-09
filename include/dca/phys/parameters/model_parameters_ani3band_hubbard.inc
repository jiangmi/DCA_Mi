// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization of model parameters for three-band Hubbard model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::ani3band_hubbard<PointGroup>>> {
public:
  ModelParameters() {}

  template <typename Concurrency>
  int getBufferSize(const Concurrency& concurrency) const;
  template <typename Concurrency>
  void pack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const;
  template <typename Concurrency>
  void unpack(const Concurrency& concurrency, char* buffer, int buffer_size, int& position);

  template <typename ReaderOrWriter>
  void readWrite(ReaderOrWriter& reader_or_writer);

  double get_t_pdx() const {
    return t_pdx_;
  }
  double get_t_pdy() const {
    return t_pdy_;
  }
  void set_t_pdx(const double t_pdx) {
    t_pdx_ = t_pdx;
  }
  void set_t_pdy(const double t_pdy) {
    t_pdy_ = t_pdy;
  }

  double get_t_pp() const {
    return t_pp_;
  }
  void set_t_pp(const double t_pp) {
    t_pp_ = t_pp;
  }

  double get_ep_d() const {
    return ep_d_;
  }
  void set_ep_d(const double ep_d) {
    ep_d_ = ep_d;
  }

  double get_ep_px() const {
    return ep_px_;
  }
  void set_ep_px(const double ep_px) {
    ep_px_ = ep_px;
  }

  double get_ep_py() const {
    return ep_py_;
  }
  void set_ep_py(const double ep_py) {
    ep_py_ = ep_py;
  }

  double get_U_dd() const {
    return U_dd_;
  }
  void set_U_dd(const double U_dd) {
    U_dd_ = U_dd;
  }

  double get_U_pp() const {
    return U_pp_;
  }
  void set_U_pp(const double U_pp) {
    U_pp_ = U_pp;
  }
  
  double get_V_pp() const {
    return V_pp_;
  }
  void set_V_pp(const double V_pp) {
    V_pp_ = V_pp;
  }
  
  double get_V_pp_prime() const {
    return V_pp_prime_;
  }
  void set_V_pp_prime(const double V_pp_prime) {
    V_pp_prime_ = V_pp_prime;
  }

  double get_V_pdx() const {
    return V_pdx_;
  }
  void set_V_pdx(const double V_pdx) {
    V_pdx_ = V_pdx;
  }
  
  double get_V_pdy() const {
    return V_pdy_;
  }
  void set_V_pdy(const double V_pdy) {
    V_pdy_ = V_pdy;
  }

private:
  double t_pdx_ = 0;
  double t_pdy_ = 0;
  double t_pp_ = 0;
  double ep_px_ = 0;
  double ep_py_ = 0;
  double ep_d_ = 0;
  double U_dd_ = 0;
  double U_pp_ = 0;
  double V_pp_ = 0;
  double V_pp_prime_ = 0;
  double V_pdx_ = 0;
  double V_pdy_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::ani3band_hubbard<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_pdx_);
  buffer_size += concurrency.get_buffer_size(t_pdy_);
  buffer_size += concurrency.get_buffer_size(t_pp_);
  buffer_size += concurrency.get_buffer_size(ep_px_);
  buffer_size += concurrency.get_buffer_size(ep_py_);
  buffer_size += concurrency.get_buffer_size(ep_d_);
  buffer_size += concurrency.get_buffer_size(U_dd_);
  buffer_size += concurrency.get_buffer_size(U_pp_);
  buffer_size += concurrency.get_buffer_size(V_pp_);
  buffer_size += concurrency.get_buffer_size(V_pp_prime_);
  buffer_size += concurrency.get_buffer_size(V_pdx_);
  buffer_size += concurrency.get_buffer_size(V_pdy_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::ani3band_hubbard<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_pdx_);
  concurrency.pack(buffer, buffer_size, position, t_pdy_);
  concurrency.pack(buffer, buffer_size, position, t_pp_);
  concurrency.pack(buffer, buffer_size, position, ep_px_);
  concurrency.pack(buffer, buffer_size, position, ep_py_);
  concurrency.pack(buffer, buffer_size, position, ep_d_);
  concurrency.pack(buffer, buffer_size, position, U_dd_);
  concurrency.pack(buffer, buffer_size, position, U_pp_);
  concurrency.pack(buffer, buffer_size, position, V_pp_);
  concurrency.pack(buffer, buffer_size, position, V_pp_prime_);
  concurrency.pack(buffer, buffer_size, position, V_pdx_);
  concurrency.pack(buffer, buffer_size, position, V_pdy_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::ani3band_hubbard<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_pdx_);
  concurrency.unpack(buffer, buffer_size, position, t_pdy_);
  concurrency.unpack(buffer, buffer_size, position, t_pp_);
  concurrency.unpack(buffer, buffer_size, position, ep_px_);
  concurrency.unpack(buffer, buffer_size, position, ep_py_);
  concurrency.unpack(buffer, buffer_size, position, ep_d_);
  concurrency.unpack(buffer, buffer_size, position, U_dd_);
  concurrency.unpack(buffer, buffer_size, position, U_pp_);
  concurrency.unpack(buffer, buffer_size, position, V_pp_);
  concurrency.unpack(buffer, buffer_size, position, V_pp_prime_);
  concurrency.unpack(buffer, buffer_size, position, V_pdx_);
  concurrency.unpack(buffer, buffer_size, position, V_pdy_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::ani3band_hubbard<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("threebands-Hubbard-model");

    try {
      reader_or_writer.execute("t_pdx", t_pdx_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_pdy", t_pdy_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_pp", t_pp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_d", ep_d_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_px", ep_px_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_py", ep_py_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U_dd", U_dd_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("U_pp", U_pp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V_pp", V_pp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V_pp_prime", V_pp_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("V_pdx", V_pdx_);
    }
    catch (const std::exception& r_e) {
    }  
    try {
      reader_or_writer.execute("V_pdy", V_pdy_);
    }
    catch (const std::exception& r_e) {
    }  

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}