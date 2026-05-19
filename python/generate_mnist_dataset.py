from enum import Enum
from io import BufferedReader
import numpy as np


class DataType(Enum):
    UnsignedByte = 0x8
    SignedByte = 0x9
    Short = 0xB
    Int = 0xC
    Float = 0xD
    Double = 0xE

    def read(self, file) -> int | float:
        match self:
            case self.UnsignedByte:
                return int.from_bytes(file.read(1))
            case self.SignedByte:
                return int.from_bytes(file.read(1), signed=True)
            case self.Short:
                return int.from_bytes(file.read(2), signed=True)
            case self.Int:
                return int.from_bytes(file.read(4), signed=True)
            case self.Float:
                return float.fromhex(file.read(4).hex())
            case self.Double:
                return float.fromhex(file.read(8).hex())
            case _:
                raise ValueError("Should be impossible")


class MNIST:
    dimensions: int
    data_type: DataType
    dimensions_sizes: list[int]

    def __init__(self, filepath: str) -> None:
        with open(filepath, "br") as file:
            self.read_metadata(file)
            self.data = np.zeros(shape=self.dimensions_sizes)
            with np.nditer(self.data, op_flags=["readwrite"]) as it:
                for x in it:
                    x[...] = self.data_type.read(file)
            # for x in range(28):
            #     print([int(y) for y in self.data[1, x, :]])

    def read_metadata(self, file: BufferedReader):
        self.read_magic_number(file)
        self.read_dimensions(file)

    def read_magic_number(self, file: BufferedReader):
        # first two bytes are always zero
        file.read(2)

        # third byte is the data type
        data_type_byte = DataType.UnsignedByte.read(file)
        if data_type_byte == 0x8:
            self.data_type = DataType.UnsignedByte
        elif data_type_byte == 0x9:
            self.data_type = DataType.SignedByte
        elif data_type_byte == 0xB:
            self.data_type = DataType.Short
        elif data_type_byte == 0xC:
            self.data_type = DataType.Int
        elif data_type_byte == 0xD:
            self.data_type = DataType.Float
        elif data_type_byte == 0xE:
            self.data_type = DataType.Double
        else:
            raise ValueError("Bad byte detected: ", data_type_byte)

        # forth byte is the amount of dimensions
        self.dimensions = int(DataType.UnsignedByte.read(file))

    def read_dimensions(self, file: BufferedReader):
        self.dimensions_sizes = [
            int(DataType.Int.read(file)) for _ in range(self.dimensions)
        ]

    def cut_edge(self):
        new_shape = (self.data.shape[0], self.data.shape[1] - 2, self.data.shape[2] - 2)
        new_data = np.zeros(shape=new_shape)
        new_data[:, :, :] = self.data[:, 1:-1, 1:-1]
        self.data = new_data

    def resample(self):
        new_shape = (
            self.data.shape[0],
            self.data.shape[1] // 2,
            self.data.shape[2] // 2,
        )
        new_data = np.zeros(shape=new_shape)
        folded_data = np.array(
            [
                self.data[:, ::2, ::2],
                self.data[:, 1::2, ::2],
                self.data[:, ::2, 1::2],
                self.data[:, 1::2, 1::2],
            ]
        )
        new_data[:, :, :] = np.average(folded_data, axis=0)
        self.data = new_data

    def show(self):
        for x in range(self.data.shape[1]):
            print([int(y) for y in self.data[1, x, :]])

    def binary(self, threshold: int = 1):
        self.data[self.data < threshold] = 0
        self.data[self.data >= threshold] = 1

    def save(self, filepath: str):
        output = []
        for x in range(self.data.shape[0]):
            numbers = self.data[x, :, :]
            string = "".join([str(int(n)) for n in np.nditer(numbers)])
            string += "\n"
            output.append(string)
        output.sort()
        with open(filepath, "w") as file:
            file.writelines(output)


if __name__ == "__main__":
    mnist = MNIST("./python/train-images-idx3-ubyte")
    # mnist.cut_edge()
    # mnist.cut_edge()
    # mnist.cut_edge()
    # mnist.resample()
    mnist.binary(threshold=100)
    mnist.show()
    mnist.save("./python/MNIST28.sorted")
