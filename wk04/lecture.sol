// Solidity Code from the lecture

contract C {
    function add(uint256 x1, uint256 y1, uint256 x2, uint256 y2) external view returns (uint256, uint256) {
        (bool ok, bytes memory result) = address(6).staticcall(abi.encode(x1, y1, x2, y2));
        require(ok, "call failed");
        return abi.decode(result, (uint256, uint256));
    }

    function multiply(uint256 x, uint256 y, uint256 s) external view returns (uint256, uint256) {
        (bool ok, bytes memory result) = address(7).staticcall(abi.encode(x, y, s));
        require(ok, "call failed");
        return abi.decode(result, (uint256, uint256));
    }
}
