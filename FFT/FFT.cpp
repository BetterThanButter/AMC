#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>

// Реализуте пропущенные методы
// в качестве Т будет использоваться std::complex<float> // <double> // <long double>
namespace FFT {
    // Возвращает корень степени degree из 1
    template <typename T>
    T GetRoot(size_t degree);

    // Выполняет преобразование фурье квадратичным алгоритмом для вектора произвольной длины
    template <typename T>
    std::vector<T> FourierTransform(const std::vector<T>& data);

    // Выполняет обратное преобразование квадратичным алгоритмом для вектора фиксированной длины
    template <typename T>
    std::vector<T> InverseFourierTransform(const std::vector<T>& data);

    // Добивает вектор в конце нулями до длины expected_length,
    // выбрасывает std::runtime_error если expected_length < data.size()
    template <typename T>
    std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length);

    // Быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastFourierTransform(const std::vector<T>& data);

    // Обратное быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastInverseFourierTransform();
} // namespace FFT


// Операции над многочленами с помощью ффт
template <typename T>
class Polynomial {
public:
    explicit Polynomial(const std::vector<T>& coefficients)
        : coefficients_(coefficients), degree_(coefficients.size()) {
    }

    // Чтобы можно было написать Polynomial<std::complex<double>> p({1, 2, 1})
    // И получить представление многочлена 1 + 2x + x^2
    Polynomial(const std::initializer_list<T>& coefficients)
        : coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size()) {
    }

    Polynomial(const Polynomial&) = default;
    Polynomial(Polynomial&&) noexcept = default;

    Polynomial& operator=(const Polynomial) = default;
    Polynomial& operator=(Polynomial&&) noexcept = default;

    size_t GetDegree() const {
        return degree_;
    }

    // Реализуйте следующие методы
    bool operator==(const Polynomial& other);

    Polynomial& operator+=(const Polynomial& other);

    Polynomial& operator-=(const Polynomial& other);

    // Возведение в степень pow с помощью комбинации FFT и индийского возведения в степень
    Polynomial& operator^=(size_t pow);

    // С Использованием быстрого преобразования фурье
    Polynomial& operator*=(const Polynomial& other);

    friend std::ostream& operator<<(std::ostream& ostr, const Polynomial&);

    // Используя предыдущее
    Polynomial operator+(const Polynomial& other);

    Polynomial operator-(const Polynomial& other);

    Polynomial operator*(const Polynomial& other);

    Polynomial operator^(size_t pow);

    // И еще один, унарный минус
    friend Polynomial operator-(const Polynomial& other);

private:
    std::vector<T> coefficients_;
    size_t degree_;
};


// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов



// Задачи, решаемые с помощью умножения многочленов
// Если вы напишете решение, работающее для произольных строк над ascii кодировкой - укажете это и
// возможно, получите небольшой бонусный балл
namespace SubstringMatching {

    // Метод принимает две строки str и pattern, возвращает индексов,
    // указывающих начала подстрок str, совпадающих с pattern
    // Можете считать, что str и pattern состоят из символов 0 и 1
    std::vector<size_t> FindSubstrings(const std::string& str, const std::string& pattern);

    // Аналогично предыдущему, но теперь pattern может содержать символ '?', на месте которого
    // можно поставить любой символ алфавита
    std::vector<size_t> FindMatches(const std::string& str, const std::string& pattern);

} // namespace SubstringMatcher



// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов
// Обратите внимание, что при проведении прямого и обратного преобразования для многочлена с
// целыми коэффициентами, новые коэффициенты могут выйти дробными и трубующими округления.