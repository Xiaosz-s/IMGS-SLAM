import lietorch_backends
import torch
import torch.nn.functional as F


class GroupOp(torch.autograd.Function):

    @classmethod
    def forward(cls, ctx, group_id, *inputs):
        ctx.group_id = group_id
        ctx.save_for_backward(*inputs)
        out = cls.forward_op(ctx.group_id, *inputs)
        return out

    @classmethod
    def backward(cls, ctx, grad):
        error_str = "Backward operation not implemented for {}".format(cls)
        assert cls.backward_op is not None, error_str
        inputs = ctx.saved_tensors
        grad = grad.contiguous()
        grad_inputs = cls.backward_op(ctx.group_id, grad, *inputs)
        return (None,) + tuple(grad_inputs)


class Exp(GroupOp):
    (forward_op, backward_op) = (
        lietorch_backends.expm,
        lietorch_backends.expm_backward,
    )


class Log(GroupOp):
    (forward_op, backward_op) = (
        lietorch_backends.logm,
        lietorch_backends.logm_backward,
    )


class Inv(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.inv, lietorch_backends.inv_backward)


class Mul(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.mul, lietorch_backends.mul_backward)


class Adj(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.adj, lietorch_backends.adj_backward)


class AdjT(GroupOp):
    (forward_op, backward_op) = (
        lietorch_backends.adjT,
        lietorch_backends.adjT_backward,
    )


class Act3(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.act, lietorch_backends.act_backward)


class Act4(GroupOp):
    (forward_op, backward_op) = (
        lietorch_backends.act4,
        lietorch_backends.act4_backward,
    )


class Jinv(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.Jinv, None)


class ToMatrix(GroupOp):
    (forward_op, backward_op) = (lietorch_backends.as_matrix, None)


class FromVec(torch.autograd.Function):

    @classmethod
    def forward(cls, ctx, group_id, *inputs):
        ctx.group_id = group_id
        ctx.save_for_backward(*inputs)
        return inputs[0]

    @classmethod
    def backward(cls, ctx, grad):
        inputs = ctx.saved_tensors
        J = lietorch_backends.projector(ctx.group_id, *inputs)
        return (
            None,
            torch.matmul(grad.unsqueeze(-2), torch.linalg.pinv(J)).squeeze(-2),
        )


class ToVec(torch.autograd.Function):

    @classmethod
    def forward(cls, ctx, group_id, *inputs):
        ctx.group_id = group_id
        ctx.save_for_backward(*inputs)
        return inputs[0]

    @classmethod
    def backward(cls, ctx, grad):
        inputs = ctx.saved_tensors
        J = lietorch_backends.projector(ctx.group_id, *inputs)
        return (None, torch.matmul(grad.unsqueeze(-2), J).squeeze(-2))
