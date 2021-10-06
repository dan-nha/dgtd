#ifndef GET_H
#define GET_H

namespace Get{

  template<typename T>
    const T left(const std::tuple<T,T> arg) {
      return std::get<0>(arg);
    }

  template<typename T>
    const T right(const std::tuple<T,T> arg) {
      return std::get<1>(arg);
    }

} // namespace Get

#endif
