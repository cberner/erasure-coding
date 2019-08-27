use crate::gf256::{Polynomial, Octet};

#[derive(Clone, Debug, PartialEq)]
pub struct BlockPolynomial {
    // largest is 0: x^2 + x + c
    // x^2 is 0 index
    // vec of coefficients. First one is all of the x^2 coefficients
    pub coefficient_arrays: Vec<Vec<u8>>
}

impl BlockPolynomial {
    pub fn new(coefficient_arrays: Vec<Vec<u8>>) -> BlockPolynomial {
        let length = coefficient_arrays[0].len();
        assert!(coefficient_arrays.iter().map(Vec::len).all(|x| x == length));
        BlockPolynomial {
            coefficient_arrays: coefficient_arrays as Vec<Vec<u8>>
        }
    }

    pub fn div(&self, divisor: &Polynomial) -> (BlockPolynomial, BlockPolynomial) {
        let mut result = self.coefficient_arrays.clone();

        for position in 0..self.coefficient_arrays[0].len() {
            for i in 0..(self.coefficient_arrays.len() - (divisor.coefficients.len() - 1)) {
                let coefficient = result[i][position].clone();
                for j in 1..divisor.coefficients.len() {
                    result[i + j][position] = (Octet::new(result[i + j][position]) + (&divisor.coefficients[j] * &Octet::new(coefficient))).byte();
                }
            }
        }
        let separator = result.len() - (divisor.coefficients.len() - 1);

        let quotient = BlockPolynomial {
            coefficient_arrays: result[..separator].to_owned()
        };
        let remainder = BlockPolynomial{
            coefficient_arrays: result[separator..].to_owned()
        };

        return (quotient, remainder);
    }
}
