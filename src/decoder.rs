use crate::gf256::{Polynomial, Octet};

pub struct Decoder {
    data_blocks: u8,
    repair_blocks: u8
}

impl Decoder {
    pub fn new(data_blocks: u8, repair_blocks: u8) -> Decoder {
        Decoder {
            data_blocks,
            repair_blocks
        }
    }

    fn calculate_syndrome(&self, data: &[Option<u8>], repair: &[u8]) -> Polynomial {
        let mut full_with_zeros: Vec<u8> = data.iter().map(|x| x.unwrap_or(0)).collect();
        full_with_zeros.extend(repair);
        let full_poly = Polynomial::new(&full_with_zeros);
        // Pad with a zero
        let mut syndrome = vec![0; self.repair_blocks as usize + 1];
        for i in 0..self.repair_blocks {
            syndrome[i as usize + 1] = full_poly.eval(&Octet::alpha(i)).byte();
        }
        return Polynomial::new(&syndrome);
    }

    fn calculate_erasure_locator(&self, data: &[Option<u8>]) -> Polynomial {
        let mut locator = Polynomial::new(&[1]);
        for i in 0..data.len() {
            if data[i].is_none() {
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                locator = locator.mul(&Polynomial::new(&[Octet::alpha(location as u8).byte(), 1]));
            }
        }
        return locator;
    }

    fn calculate_erasure_evaluator(&self, syndrome: &Polynomial, erasure_locator: &Polynomial) -> Polynomial {
        let mut poly = Polynomial::new(&vec![0; (erasure_locator.coefficients.len() + 1) as usize]);
        poly.coefficients[0] = Octet::one();

        let (_, remainder) = syndrome.mul(erasure_locator).div(&poly);
        return remainder
    }

    // Forney algorithm
    fn correct_erasures(&self, syndrome: &Polynomial, data: &[Option<u8>], repair: &[u8]) -> Polynomial {
        let locator = self.calculate_erasure_locator(data);
        let mut syndrome = syndrome.clone();
        syndrome.coefficients.reverse();
        let mut erasure_evaluator = self.calculate_erasure_evaluator(&syndrome, &locator);
        erasure_evaluator.coefficients.reverse();

        let mut X = vec![];
        let mut full_message_erasure_positions = vec![];
        for i in 0..data.len() {
            if data[i].is_none() {
                full_message_erasure_positions.push(i);
                let location = self.data_blocks + self.repair_blocks - 1 - i as u8;
                let l = 255 - location as i32;
                X.push(Octet::new(2).pow(-l));
            }
        }

        let mut E = vec![Octet::zero(); (self.data_blocks + self.repair_blocks) as usize];
        let X_length = X.len();
        for (i, Xi) in X.iter().enumerate() {
            let Xi_inv = &Octet::one() / Xi;

            let mut prime_tmp = vec![];
            for j in 0..X_length {
                if j != i {
                    prime_tmp.push(Octet::one() - &Xi_inv * &X[j]);
                }
            }

            let mut locator_prime = Octet::one();
            for value in prime_tmp.iter() {
                locator_prime = &locator_prime * value;
            }

            erasure_evaluator.coefficients.reverse();
            let y = erasure_evaluator.eval(&Xi_inv);
            erasure_evaluator.coefficients.reverse();
            // TODO: uh, this pow() seems to do nothing...
            let y = Xi.pow(1) * y;

            assert_ne!(locator_prime, Octet::zero());
            let magnitude = y / locator_prime;

            E[full_message_erasure_positions[i]] = magnitude;
        }

        let mut full_with_zeros: Vec<u8> = data.iter().map(|x| x.unwrap_or(0)).collect();
        full_with_zeros.extend(repair);
        let full_poly = Polynomial::new(&full_with_zeros);
        let E_bytes: Vec<u8> = E.iter().map(Octet::byte).collect();
        return full_poly.add(&Polynomial::new(&E_bytes));
    }

    // TODO: support missing repair blocks
    // Returns tuple of (data blocks, repair blocks)
    pub fn decode(&self, data: &[Option<u8>], repair: &[u8]) -> Vec<u8> {
        // TODO: support blocks with more than 1 byte
        assert_eq!(data.len(), self.data_blocks as usize);
        assert_eq!(repair.len(), self.repair_blocks as usize);

        let syndrome = self.calculate_syndrome(data, repair);
        let poly = self.correct_erasures(&syndrome, data, repair);
        let mut repaired = poly.into_coefficients();
        repaired.truncate(self.data_blocks as usize);

        return repaired;
    }
}

#[cfg(test)]
mod tests {
    use crate::Decoder;
    use crate::gf256::Polynomial;

    #[test]
    fn erasure_locator() {
        let data = vec![None, Some(1), Some(2), None, None];
        assert_eq!(Polynomial::new(&[19, 28, 152, 1]),
                   Decoder::new(5, 3).calculate_erasure_locator(&data));
    }
}
