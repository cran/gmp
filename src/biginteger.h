
#include <gmp.h>
#include <R.h>
#include <Rdefines.h>

#include <string>

/**
 * A big integer. Actually a wrapper for mpz_t to work with plus 
 * some special stuff.
 *
 * The biginteger special state "NA" means, no value is assigned. 
 * This does not mean, the internal state is not constructed, but 
 * the value explicit is "not available".
 */




class biginteger
{
    /**
     * The actual integer value.
     */
    mpz_t value;

    /**
     * True, if the value is "NA".
     */
    bool na;
    
public:
    /**
     * Construct a "NA" biginteger.
     */
    biginteger() {mpz_init(value); setValue();}

    /**
     * Construct a biginteger from a raw expression. 
     */
    biginteger(void* raw);

    /**
     * Create a biginteger from a value. Remember to free the 
     * parameter's mpz_t if you allocated them by yourself - 
     * biginteger will copy the value.
     */
    biginteger(const mpz_t& value_) : na(false) {mpz_init_set(value, value_);}
    
    /**
     * Construct a biginteger from a long value.
     */
    biginteger(int value_) : na(false) {mpz_init_set_ui(value, value_);}
    
    /**
     * Construct a biginteger from a double value.
     */
    biginteger(double value_) : na(false) {mpz_init_set_d(value, value_);  }

    /**
     * Construct a biginteger from a string value.
     */
    biginteger(const std::string& value_) : na(false) {mpz_init_set_str(value, value_.c_str(), 0);}
    
    /**
     *  Copy constructor (mpz_t aren't standard-copyable)
     */
    biginteger(const biginteger& rhs) : na(rhs.na) {mpz_init_set(value, rhs.value);}

    /**
     *  Assignment operator.
     */
    biginteger& operator=(const biginteger& rhs) {mpz_set(value, rhs.value); na = rhs.na;}

    /**
     * Free the owned mpz_t structs
     */
    ~biginteger() {mpz_clear(value);}


    /**
     * Set the biginteger to state "NA".
     */
    void setValue() {mpz_set_ui(value, 0); na = true;}

    /**
     * Set the biginteger to a specific value.
     */
    void setValue(mpz_t value_) {mpz_set(value, value_); na = false;}
    void setValue(int value_) {mpz_set_ui(value, value_); na = false;}

    /**
     * For const-purposes, return the value. Remember, that the return value
     * only lives as long as this class live, so do not call getValueTemp on
     * temporary objects.
     */
    const mpz_t& getValueTemp() const {return value;}

    /**
     * Return true, if the value is NA.
     */
    bool isNA() const {return na;}
    
    /**
     * Return true, if the value is 0.
     */
    int sgn() const {return mpz_sgn(value);}
    
    /**
     *  Convert the biginteger into a standard string.
     */
    std::string str() const;

    /**
     * Convert the biginteger into a long value (cut off the msb's if it don't
     * fit).
     */
    long as_long() const {return mpz_get_ui(value);}

    /**
     * \brief Convert the biginteger into a double value 
     * (you may loose precision)
     */
    double as_double() const {return mpz_get_d(value);}

    /**
     * Convert the biginteger to a raw memory block. Obtain the size needed
     * from biginteger_raw_size() first and make sure, the buffer provided is
     * large enough to hold the data.
     *
     * Also remember, that the modulus is not saved this way. To obtain a
     * modulus raw byte use get_modulus().as_raw(void*).
     * 
     * @return number of bytes used (same as raw_size())
     */
    int as_raw(void* raw) const;

    /**
     * Return the number of bytes needed to store this biginteger in a
     * continous memory block.
     */
    size_t raw_size() const;

    /**
     * Swap values with the argument
     */
    void swap(biginteger& other);


    /** 
     * Test prime numbers
     */
    int isprime(int reps){return  mpz_probab_prime_p(value,reps);}

};


/**
 * Represents two biginteger: a value and a modulus. These both are used
 * to operate arithmetic functions on it. If the modulus is NA, no modulus
 * to the operation result is applied. If the value is NA, the result is always NA.
 */
class bigmod {
public:
    biginteger value;
    biginteger modulus;
    
    /**
     * Construct a bigmod from a biginteger, optional with a given modulus
     */
    bigmod(const biginteger& value_ = biginteger(), 
	    const biginteger& modulus_ = biginteger()) : value(value_),modulus(modulus_) {}
    
    /**
     * Return as a human readible string
     */
    std::string str() const;
};




/**
 * Add two bigmods together.
 *
 * If only one has a modulus set, the result will have this
 * modulus. If both bigmods disagree with the modulus, the result will not have
 * a modulus set. If none modulus for either bigmod is set, the result will not
 * have a modulus as well.
 */
bigmod operator+(const bigmod& rhs, const bigmod& lhs);

/**
 * Subtract two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator-(const bigmod& rhs, const bigmod& lhs);

/**
 * Multiply two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator*(const bigmod& rhs, const bigmod& lhs);

/**
 * Divide two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator/(const bigmod& rhs, const bigmod& lhs);

/**
 * Calculate the modulus (remainder) of two bigmods.
 *
 * The resulting bigmod will have set the intern modulus to 
 * the value of lhs, no matter what rhs.modulus or lhs.modulus
 * was before, except if rhs and lhs has both no modulus set,
 * in which case the resulting modulus will be unset too.
 */
bigmod operator%(const bigmod& rhs, const bigmod& lhs);

/**
 * Return the power of "exp" to the base of "base" (return = base^exp).
 *
 * If both moduli are unset or unequal, this may EAT your memory alive, 
 * since then the infinite "pow" is used instead of the modulus "powm". 
 * You  may not try to pow a value this way with an exponent that does
 * not fit into a long value.
 *
 * For other modulus description, see operator+(bigmod, bigmod)
 */
bigmod pow(const bigmod& base, const bigmod& exp);

/**
 * Return the modulo inverse to x mod m. (return = x^-1 % m)
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod inv(const bigmod& x, const bigmod& m);

/**
 * Return a bigmod with value (x % m) and the intern modulus set to m.
 * Intern modulus settings of x and m are ignored.
 *
 * Do not confuse this with operator%(bigmod, bigmod).
 */
bigmod set_modulus(const bigmod& x, const bigmod& m);

/**
 * Return the greatest common divisor of both parameters
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod gcd(const bigmod& rhs, const bigmod& lhs);

/**
 * Return the least common multiply of both parameter.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod lcm(const bigmod& rhs, const bigmod& lhs);


