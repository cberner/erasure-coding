use crate::gf256::{Polynomial, fused_addassign_mul_scalar};

fn get_both_indices<T>(vector: &mut Vec<T>, i: usize, j: usize) -> (&mut T, &mut T) {
    debug_assert_ne!(i, j);
    debug_assert!(i < vector.len());
    debug_assert!(j < vector.len());
    if i < j {
        let (first, last) = vector.split_at_mut(j);
        return (&mut first[i], &mut last[0]);
    } else {
        let (first, last) = vector.split_at_mut(i);
        return (&mut last[0], &mut first[j]);
    }
}

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

        for i in 0..(self.coefficient_arrays.len() - (divisor.coefficients.len() - 1)) {
            for j in 1..divisor.coefficients.len() {
                let both: (&mut Vec<u8>, &mut Vec<u8>) = get_both_indices(&mut result, i + j, i);
                let (dest, src) = both;
                fused_addassign_mul_scalar(dest, src, &divisor.coefficients[j]);
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
