#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>
#include <cmath>

const double pi = 3.1415;

// Реализуте пропущенные методы
// в качестве Т будет использоваться std::complex<float> // <double> // <long double>
namespace FFT {
    // Возвращает корень степени degree из 1
    template <typename T>
    T GetRoot(size_t degree) {
        return T(cos(2*pi / degree), sin(2*pi / degree));
    };

    // Выполняет преобразование фурье квадратичным алгоритмом для вектора произвольной длины
    template <typename T>
    std::vector<T> FourierTransform(const std::vector<T>& data){

        std::vector<T> result;
        size_t n(data.size());
        T value(0);
        T root_n(GetRoot<T>(n));

        for (unsigned int i = 0; i < n; i++) {

            for (unsigned int j = 0; j < n; j++) {
                value += pow(root_n, i*j) * data[j];
            }
            result.push_back(value);
            value = {0,0};
        }

        return result;
    };

    // Выполняет обратное преобразование квадратичным алгоритмом для вектора фиксированной длины
    template <typename T>
    std::vector<T> InverseFourierTransform(const std::vector <T> &data) {

        size_t n(data.size());
        std::vector<T> data_returned(n);

        T value(0);
        T root_n(GetRoot<T>(n));

        T one(1);

        for (size_t i = 0; i < n; i++) {

            for (size_t j = 0; j < n; j++) {
                //bullshit code
                value += one / (pow(root_n, (i * j))) * data[j];
            }
            //normed
            data_returned[i] = (value / T(n));
            value = {0,0};
        }

        return data_returned;
    }
    // Добивает вектор в конце нулями до длины expected_length,
    // выбрасывает std::runtime_error если expected_length < data.size()
    template <typename T>
    std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length){

        std::vector<T> result(data);
        for (size_t i = data.size(); i < expected_length; i++)
            result.push_back(T(0));
        return result;
    };

    // Быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastFourierTransform(const std::vector<T>& data){

        size_t n = data.size();
        std::vector<T> result(n);
        if (n == 1){
           return data;
        }

        std::vector<T> even;
        std::vector<T> odd;

        for (size_t i = 0; i < n ; i++) {
            if(i % 2 == 0) {
                even.push_back(data[i]);
            }
            else {
                odd.push_back(data[i]);
            }
        }

        even = FastFourierTransform(even);
        odd = FastFourierTransform(odd);

        T bias = 1;
        T root_n = GetRoot<T>(n);
        for (size_t i = 0; i < n/2; i++) {
            result[i] = even[i] + bias * odd[i];
            result[i+n/2] = even[i] - bias * odd[i];

            bias *= root_n;
        }

        return result;
    }

    // Обратное быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastInverseFourierTransform(const std::vector<T>& data){
        size_t n = data.size();
        std::vector<T> result(n);
        if (n == 1){
            return data;
        }

        std::vector<T> even;
        std::vector<T> odd;

        for (size_t i = 0; i < n ; i++) {
            if(i % 2 == 0) {
                even.push_back(data[i]);
            }
            else {
                odd.push_back(data[i]);
            }
        }


        even = FastInverseFourierTransform(even);
        odd = FastInverseFourierTransform(odd);

        T bias = 1;
        T root_n = GetRoot<T>(n);
        for (size_t i = 0; i < n/2; i++) {
            //pow in negative
            //Так как размер 2^n
          //  std::cout << T(0.5) * T(2) << std::endl;
            result[i] = T(0.5)*(even[i] + pow(bias, -1) * odd[i]);
            result[i+n/2] = T(0.5)*(even[i] - pow(bias, -1) * odd[i]);

            bias *= root_n;
        }

        return result;
    }
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

    //конструктор по умолчанию
    Polynomial(const Polynomial&) = default;
    Polynomial(Polynomial&&) noexcept = default;

    //Polynomial& operator=(const Polynomial) = default;
    Polynomial& operator=(Polynomial&&) noexcept = default;

    size_t GetDegree() const {
        return degree_;
    }

    // Реализуйте следующие методы
    bool operator==(const Polynomial& other) {

        bool is_equal = false;
        if (this->coefficients_ == other.coefficients_) {
            is_equal = true;
        }
        else {
            is_equal = false;
        }
        return is_equal;
    };

    void print() {
        std::cout << "Printing polynom: " << std::endl;
        for (std::vector<std::complex<double>>::iterator i = this->begin(); i != this->end(); ++i)
            std::cout << *i << ' ';
        std::cout << std::endl;
    }

    Polynomial& operator+=(const Polynomial& other) {

        if (this->degree_ < other.degree_)
            this->coefficients_.resize(other.degree_);

        for (size_t i = 0; i < other.degree_; i++) {
            this->coefficients_[i] += other.coefficients_[i];
        }
    };

    Polynomial& operator-=(const Polynomial& other) {

        if (this->degree_ < other.degree_)
            this->coefficients_.resize(other.degree_);

        for (size_t i = 0; i < other.degree_; i++) {
            this->coefficients_[i] -= other.coefficients_[i];
        }
    };

    // Возведение в степень pow с помощью комбинации FFT и индийского возведения в степень
    Polynomial& operator^=(size_t pow){

        if (pow == 0){
            std::vector<T> v(1);
            Polynomial *e = new Polynomial(v);
            return *e;
        }
        if ( pow == 1) {
            return this;
        }

        Polynomial start_polynom = *this;

        while (pow > 0){
           if (pow % 2 == 1){
                *this *= start_polynom;
           }
            start_polynom *= start_polynom;
           pow /= 2;
        }
        return *this;
    };

    // С Использованием быстрого преобразования фурье
    Polynomial& operator*=(const Polynomial& other){
        *this = *this * other;
        return *this;
    };

    friend std::ostream& operator<<(std::ostream& ostr, const Polynomial& p) {
        //std::vector<std::complex <double >>:: iterator it;
        size_t degree = p.GetDegree();

        size_t i;
        for(i = 0; i < degree - 1; i++ )    {
            ostr << p.coefficients_[i] << "X^" << i  << " +";
        }
        ostr << p.coefficients_[i] << "X^" << i  << "\n";
        ostr << "kek";
        ostr << "\n";
        return ostr;
    };

    // Используя предыдущее
    Polynomial operator+(const Polynomial& other) {
        Polynomial answer (*this);
        answer += other;
        return answer;
    };

    Polynomial operator-(const Polynomial& other) {
        Polynomial answer (*this);
        answer -= other;
        return answer;
    };

    Polynomial operator*(const Polynomial& other) {

        std::vector<T> other_new = other.coefficients_;
        std::vector<T> this_new = this->coefficients_;

        size_t degr_sum = this->degree_ + other.degree_;
        size_t extended_size = (size_t) pow(2, std::ceil(log(degr_sum + 1)/log(2)));

        //FTT for polynoms
        std::vector<T> other_FTT = FFT::FastFourierTransform(FFT::AddPadding(other_new, extended_size));
        std::vector<T> this_FTT = FFT::FastFourierTransform(FFT::AddPadding(this_new, extended_size));

        //result
        std::vector<T> result;
        for (int i = 0; i < extended_size; i++){
            result.push_back(other_FTT[i] * this_FTT[i]);
        }

        return Polynomial(FFT::FastInverseFourierTransform(result));
    };

    Polynomial operator^(size_t pow);

    // И еще один, унарный минус
//    friend Polynomial operator-(const Polynomial& other);

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
int main() {

    const std::vector <int> foo = {3,4,5};

// Чтобы можно было написать Polynomial<std::complex<double>> p({1, 2, 1})

    Polynomial <std::complex <double >> p1({{1,0},{2,0}, {3,0}});
    Polynomial <std::complex <double >> p2({{1,0},{2,0},{3,0}});

    std::vector<std::complex<double>> a ;
    a.push_back({1,0});
    a.push_back({2,0});
    a.push_back({3,0});
    a.push_back({4,0});

    Polynomial <std::complex <double >> p3 = p1*p2;
    std::cout << p3;

    std::cout << "polynim:" << std::endl;
    for (std::vector<std::complex<double>>::iterator i = a.begin(); i != a.end(); ++i)
        std::cout << *i << ' ';

    std::vector<std::complex<double>> c3 = FFT::FastFourierTransform(a);


    std::cout << std::endl<<  " After FFFT: " << std::endl;
    for (std::vector<std::complex<double>>::const_iterator i = c3.begin(); i != c3.end(); ++i)
        std::cout << *i << ' ';
    std::cout << std::endl<<  " After IFFFT: " << std::endl;
    std::vector<std::complex<double>> b = FFT::FastInverseFourierTransform(c3);
    for (int i = 0; i < b.size(); ++i)
        std::cout << b[i] << ' ';
    return 0;

}