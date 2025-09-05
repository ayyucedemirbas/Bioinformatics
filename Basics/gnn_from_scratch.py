import numpy as np
np.random.seed(0)
def relu(x):
    return np.maximum(0, x)

class MeanAggLayerNP:
    def __init__(self, in_dim, out_dim, activation=relu):
        self.W = np.random.randn(in_dim, out_dim) * np.sqrt(2.0 / in_dim)
        self.B = np.random.randn(in_dim, in_dim) * np.sqrt(2.0 / in_dim)
        self.activation = activation
        # backprop
        self.dW = np.zeros_like(self.W)
        self.dB = np.zeros_like(self.B)

    def forward(self, H, A):
        deg = np.maximum(A.sum(axis=1, keepdims=True), 1e-10)
        neighbor_mean = (A @ H) / deg
        self.H = H  #backprop
        self.A = A
        self.neighbor_mean = neighbor_mean

        combined = neighbor_mean + (H @ self.B.T)
        self.combined = combined
        Z = combined @ self.W
        self.Z = Z
        return self.activation(Z) if self.activation else Z

    def backward(self, grad_out):
        dZ = grad_out * (self.Z > 0).astype(float)  # ReLU grad
        self.dW = self.combined.T @ dZ
        d_combined = dZ @ self.W.T

        dH_self = d_combined @ self.B
        deg = np.maximum(self.A.sum(axis=1, keepdims=True), 1e-10)
        dH_neigh = (self.A.T @ (d_combined / deg))  # approximate

        return dH_self + dH_neigh

class SimpleGCNNP:
    def __init__(self, in_dim, hidden_dim, out_dim, num_layers=2, lr=0.01):
        dims = [in_dim] + [hidden_dim]*(num_layers-1) + [out_dim]
        self.layers = [MeanAggLayerNP(dims[i], dims[i+1], activation=(relu if i < len(dims)-2 else None))
                       for i in range(len(dims)-1)]
        self.lr = lr

    def forward(self, X, A):
        H = X
        for layer in self.layers:
            H = layer.forward(H, A)
        return H

    def backward(self, grad):
        for layer in reversed(self.layers):
            grad = layer.backward(grad)

    def step(self):
        for layer in self.layers:
            layer.W -= self.lr * layer.dW
            layer.B -= self.lr * layer.dB

    
A = np.array([
    [0,1,1,0,0,0],
    [1,0,1,0,0,0],
    [1,1,0,1,0,0],
    [0,0,1,0,1,1],
    [0,0,0,1,0,1],
    [0,0,0,1,1,0],
], dtype=float)

X = np.random.randn(6, 8)
y = np.array([0,0,0,1,1,1])

model = SimpleGCNNP(8, 16, 2, num_layers=2, lr=0.01)
for epoch in range(300):
    logits = model.forward(X, A)  # [6,2]
    exp = np.exp(logits - logits.max(axis=1, keepdims=True))
    probs = exp / exp.sum(axis=1, keepdims=True)

    loss = -np.mean(np.log(probs[np.arange(6), y] + 1e-9))
    if epoch % 50 == 0:
        preds = probs.argmax(axis=1)
        acc = np.mean(preds == y)
        print(f"epoch {epoch:03d} | loss {loss:.4f} | acc {acc:.3f}")

    grad_logits = probs
    grad_logits[np.arange(6), y] -= 1
    grad_logits /= 6
    model.backward(grad_logits)
    model.step()

    print(probs)
