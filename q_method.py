import numpy as np

def q_method(b_vectors, r_vectors, weights=None):
    """
    Implementation of the q-method for Wahba's problem.

    Parameters:
    - b_vectors: ndarray of shape (N, 3)
        Array of N vectors measured in the body frame.
    - r_vectors: ndarray of shape (N, 3)
        Array of N reference vectors in the inertial frame.
    - weights: ndarray of shape (N,), optional
        Array of weights for each vector pair. Default is equal weighting.

    Returns:
    - q_opt: ndarray of shape (4,)
        Optimal quaternion as [q0, q1, q2, q3] (scalar-first convention).
    - R_opt: ndarray of shape (3, 3)
        Optimal rotation matrix.
    """
    # Ensure input vectors have the correct shapes
    b_vectors = np.asarray(b_vectors)
    r_vectors = np.asarray(r_vectors)
    if weights is None:
        weights = np.ones(b_vectors.shape[0])
    else:
        weights = np.asarray(weights)

    # Validate dimensions
    assert b_vectors.shape == r_vectors.shape, "b_vectors and r_vectors must have the same shape."
    assert b_vectors.shape[1] == 3, "Each vector must have 3 components."
    assert weights.shape[0] == b_vectors.shape[0], "weights must have the same length as the number of vectors."

    # Step 1: Compute the attitude profile matrix K
    K = np.zeros((3, 3))
    for i in range(b_vectors.shape[0]):
        K += weights[i] * np.outer(b_vectors[i], r_vectors[i])

    # Step 2: Construct the 4x4 symmetric matrix Q
    trace_K = np.trace(K)
    Q = np.zeros((4, 4))
    Q[0, 0] = trace_K
    Q[1:, 1:] = K + K.T - np.eye(3) * trace_K
    Q[0, 1:] = Q[1:, 0] = [K[1, 2] - K[2, 1], K[2, 0] - K[0, 2], K[0, 1] - K[1, 0]]

    # Step 3: Compute the largest eigenvalue and corresponding eigenvector of Q
    eigenvalues, eigenvectors = np.linalg.eigh(Q)
    max_index = np.argmax(eigenvalues)
    q_opt = eigenvectors[:, max_index]  # Optimal quaternion

    # Step 4: Convert quaternion to rotation matrix (optional)
    q0, q1, q2, q3 = q_opt
    R_opt = np.array([
        [q0**2 + q1**2 - q2**2 - q3**2, 2 * (q1*q2 - q0*q3), 2 * (q1*q3 + q0*q2)],
        [2 * (q1*q2 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2 * (q2*q3 - q0*q1)],
        [2 * (q1*q3 - q0*q2), 2 * (q2*q3 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]
    ])

    return q_opt, R_opt

# Example usage
if __name__ == "__main__":
    # Define body-frame vectors
    body_vectors = np.array([
        [0.2673, 0.5345, 0.8018],
        [0.9636, 0.1483, 0.2224],
        [0.2113, 0.7887, 0.5774]
    ])

    inertial_vectors = np.array([
        [0.7071, 0.7071, 0.0],
        [0.8660, 0.5, 0.0],
        [0.0, 0.7071, 0.7071]
    ])

    # Solve for the optimal quaternion and rotation matrix
    q_opt, R_opt = q_method(body_vectors, inertial_vectors)

    # Print results
    print("Optimal Quaternion (q0, q1, q2, q3):", q_opt)
    print("Optimal Rotation Matrix:")
    print(R_opt)
