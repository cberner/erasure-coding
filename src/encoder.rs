use crate::gf256::Polynomial;
use crate::base::Block;
use crate::block_polynomial::BlockPolynomial;

pub struct Encoder {
    data_blocks: u8,
    repair_blocks: u8
}

impl Encoder {
    pub fn new(data_blocks: u8, repair_blocks: u8) -> Encoder {
        Encoder {
            data_blocks,
            repair_blocks
        }
    }

    // Returns tuple of (data blocks, repair blocks)
    pub fn encode(&self, data: &[u8]) -> (Vec<Block>, Vec<Block>) {
        assert_eq!(data.len() % self.data_blocks as usize, 0);
        let block_length = data.len() / self.data_blocks as usize;
        // Data striped across blocks: [0, block_length) in the first block,
        // then [block_length, 2 * block_length) in the second
        let mut data_blocks = vec![];
        for i in 0..self.data_blocks as usize {
            data_blocks.push(data[i * block_length..(i + 1)*block_length].to_vec())
        }

        let mut encode_blocks = data_blocks.clone();
        let repair_blocks = vec![vec![0; block_length]; self.repair_blocks as usize];
        encode_blocks.extend(repair_blocks);
        let generator_polynomial = Polynomial::create_generator_polynomial(self.repair_blocks);

        let block_poly = BlockPolynomial::new(encode_blocks);
        let (_, repair_poly) = block_poly.div(&generator_polynomial);

        return (data_blocks, repair_poly.coefficient_arrays);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn encode() {
        let data = vec![1, 2, 3];
        let expected_data = vec![vec![1], vec![2], vec![3]];
        let expected_repair = vec![vec![4], vec![4]];

        let encoder = super::Encoder::new(data.len() as u8, expected_repair.len() as u8);
        let (encoded_data, repair) = encoder.encode(&data);

        assert_eq!(expected_data, encoded_data);
        assert_eq!(expected_repair, repair);
    }
}
