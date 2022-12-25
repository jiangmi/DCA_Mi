// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific
// publications.
//
// Author: Mi Jiang
//
// Template specialization of model parameters for d-p-s model.

template <typename PointGroup>
class ModelParameters<models::TightBindingModel<models::dpsNimodel<PointGroup>>> {
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

  double get_t_pd() const {
    return t_pd_;
  }
  void set_t_pd(const double t_pd) {
    t_pd_ = t_pd;
  }

  double get_t_pp() const {
    return t_pp_;
  }
  void set_t_pp(const double t_pp) {
    t_pp_ = t_pp;
  }

  double get_ts() const {
    return ts_;
  }
  void set_ts(const double ts) {
    ts_ = ts;
  }

  double get_ts_prime() const {
    return ts_prime_;
  }
  void set_ts_prime(const double ts_prime) {
    ts_prime_ = ts_prime;
  }

  double get_t_perp() const {
    return t_perp_;
  }
  void set_t_perp(const double t_perp) {
    t_perp_ = t_perp;
  }

  double get_t_perp_prime() const {
    return t_perp_prime_;
  }
  void set_t_perp_prime(const double t_perp_prime) {
    t_perp_prime_ = t_perp_prime;
  }

  double get_ep_d() const {
    return ep_d_;
  }
  void set_ep_d(const double ep_d) {
    ep_d_ = ep_d;
  }

  double get_ep_p() const {
    return ep_p_;
  }
  void set_ep_p(const double ep_p) {
    ep_p_ = ep_p;
  }

  double get_es() const {
    return es_;
  }
  void set_es(const double es) {
    es_ = es;
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

private:
  double t_pd_ = 0;
  double t_pp_ = 0;
  double ts_ = 0;
  double ts_prime_ = 0;
  double t_perp_ = 0;
  double t_perp_prime_ = 0;
  double ep_p_ = 0;
  double ep_d_ = 0;
  double es_ = 0;
  double U_dd_ = 0;
  double U_pp_ = 0;
};

template <typename PointGroup>
template <typename Concurrency>
int ModelParameters<models::TightBindingModel<models::dpsNimodel<PointGroup>>>::getBufferSize(
    const Concurrency& concurrency) const {
  int buffer_size = 0;

  buffer_size += concurrency.get_buffer_size(t_pd_);
  buffer_size += concurrency.get_buffer_size(t_pp_);
  buffer_size += concurrency.get_buffer_size(ts_);
  buffer_size += concurrency.get_buffer_size(ts_prime_);
  buffer_size += concurrency.get_buffer_size(t_perp_);
  buffer_size += concurrency.get_buffer_size(t_perp_prime_);
  buffer_size += concurrency.get_buffer_size(ep_p_);
  buffer_size += concurrency.get_buffer_size(ep_d_);
  buffer_size += concurrency.get_buffer_size(es_);
  buffer_size += concurrency.get_buffer_size(U_dd_);
  buffer_size += concurrency.get_buffer_size(U_pp_);

  return buffer_size;
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::dpsNimodel<PointGroup>>>::pack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, t_pd_);
  concurrency.pack(buffer, buffer_size, position, t_pp_);
  concurrency.pack(buffer, buffer_size, position, ts_);
  concurrency.pack(buffer, buffer_size, position, ts_prime_);
  concurrency.pack(buffer, buffer_size, position, t_perp_);
  concurrency.pack(buffer, buffer_size, position, t_perp_prime_);
  concurrency.pack(buffer, buffer_size, position, ep_p_);
  concurrency.pack(buffer, buffer_size, position, ep_d_);
  concurrency.pack(buffer, buffer_size, position, es_);
  concurrency.pack(buffer, buffer_size, position, U_dd_);
  concurrency.pack(buffer, buffer_size, position, U_pp_);
}

template <typename PointGroup>
template <typename Concurrency>
void ModelParameters<models::TightBindingModel<models::dpsNimodel<PointGroup>>>::unpack(
    const Concurrency& concurrency, char* buffer, int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, t_pd_);
  concurrency.unpack(buffer, buffer_size, position, t_pp_);
  concurrency.unpack(buffer, buffer_size, position, ts_);
  concurrency.unpack(buffer, buffer_size, position, ts_prime_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_);
  concurrency.unpack(buffer, buffer_size, position, t_perp_prime_);
  concurrency.unpack(buffer, buffer_size, position, ep_p_);
  concurrency.unpack(buffer, buffer_size, position, ep_d_);
  concurrency.unpack(buffer, buffer_size, position, es_);
  concurrency.unpack(buffer, buffer_size, position, U_dd_);
  concurrency.unpack(buffer, buffer_size, position, U_pp_);
}

template <typename PointGroup>
template <typename ReaderOrWriter>
void ModelParameters<models::TightBindingModel<models::dpsNimodel<PointGroup>>>::readWrite(
    ReaderOrWriter& reader_or_writer) {
  try {
    reader_or_writer.open_group("dps-IL-Nickelate-model");

    try {
      reader_or_writer.execute("t_pd", t_pd_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_pp", t_pp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ts", ts_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ts_prime", ts_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_perp", t_perp_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("t_perp_prime", t_perp_prime_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_d", ep_d_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("ep_p", ep_p_);
    }
    catch (const std::exception& r_e) {
    }
    try {
      reader_or_writer.execute("es", es_);
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

    reader_or_writer.close_group();
  }
  catch (const std::exception& r_e) {
  }
}
